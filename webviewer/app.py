
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib.parse import urljoin, urlparse

import psycopg2
from psycopg2.extras import RealDictCursor
from flask import (
    Flask,
    abort,
    g,
    jsonify,
    redirect,
    render_template,
    request,
    send_from_directory,
    session,
    url_for,
)

BASE_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = (BASE_DIR.parent / "ver6_all_analysis" / "png").resolve()
DEFAULT_PAGE_SIZE = 25
MAX_PAGE_SIZE = 100
ACCESS_PASSWORD = os.environ.get("IDRCC_PASSWORD", "ShimoLAB0501")
SESSION_KEY = os.environ.get("IDRCC_SECRET_KEY", "replace-this-secret")
SUPABASE_DB_URL = os.environ.get("SUPABASE_DB_URL")
if not SUPABASE_DB_URL:
    raise RuntimeError("SUPABASE_DB_URL environment variable is required.")

PROTEINS_VER6 = "proteins_ver6"
PROTEINS_VER9 = "proteins_ver9"
PROTEINS_VER10 = "proteins_ver10"
PPI_EDGES = "ppi_edges"
IDR_SEGMENTS_VER9 = "idr_segments_ver9"

SEARCH_MODE_COLUMN_MAP: Dict[str, List[str]] = {
    "all": ["uniprot_id", "gene_name", "protein_name", "subcellular_location"],
    "uniprot": ["uniprot_id"],
    "gene": ["gene_name"],
    "protein": ["protein_name"],
    "location": ["subcellular_location"],
}
SEARCH_MODE_OPTIONS: List[Tuple[str, str]] = [
    ("all", "All fields"),
    ("uniprot", "UniProt ID"),
    ("gene", "Gene Name"),
    ("protein", "Protein Name"),
    ("location", "Subcellular Location"),
]

app = Flask(__name__)
app.secret_key = SESSION_KEY

PROTEIN_DISPLAY_COLUMNS: List[Tuple[str, str]] = [
    ("uniprot_id", "UniProt_ID"),
    ("gene_name", "Gene_Name"),
    ("protein_name", "Protein_Name"),
    ("sequence_length", "Sequence_Length"),
    ("sequence", "Sequence"),
    ("status", "Status"),
    ("is_selenoprotein", "Is_Selenoprotein"),
    ("processing_status", "Processing_Status"),
    ("unknown_aa_count", "Unknown_AA_Count"),
    ("unknown_aa_percentage", "Unknown_AA_Percentage"),
    ("has_idr", "Has_IDR"),
    ("num_idrs", "Num_IDRs"),
    ("idr_boundaries", "IDR_Boundaries"),
    ("idr_residues", "IDR_Residues"),
    ("idr_percentage", "IDR_Percentage"),
    ("mean_disorder_score", "Mean_Disorder_Score"),
    ("max_disorder_score", "Max_Disorder_Score"),
    ("num_cc_domains", "Num_CC_Domains"),
    ("total_cc_length", "Total_CC_Length"),
    ("cc_percentage", "CC_Percentage"),
    ("longest_cc_domain_length", "Longest_CC_Domain_Length"),
    ("mean_cc_domain_length", "Mean_CC_Domain_Length"),
    ("cc_mean_score", "CC_Mean_Score"),
    ("cc_max_score", "CC_Max_Score"),
    ("subcellular_location", "Subcellular_Location"),
    ("domain_information", "Domain_Information"),
    ("interactors_80pct", "Interactors_80pct"),
    ("diseases_80pct", "Diseases_80pct"),
    ("processes_80pct", "Processes_80pct"),
    ("functions_80pct", "Functions_80pct"),
    ("interactors_90pct", "Interactors_90pct"),
    ("diseases_90pct", "Diseases_90pct"),
    ("processes_90pct", "Processes_90pct"),
    ("functions_90pct", "Functions_90pct"),
]

IDR_DISPLAY_COLUMNS: List[Tuple[str, str]] = PROTEIN_DISPLAY_COLUMNS + [
    ("idr_number", "IDR_Number"),
    ("idr_start", "IDR_Start"),
    ("idr_end", "IDR_End"),
    ("idr_length", "IDR_Length"),
    ("idr_sequence", "IDR_Sequence"),
    ("idr_percentage_in_protein", "IDR_Percentage_in_Protein"),
    ("cluster_k30", "cluster_k30"),
    ("d_min", "d_min"),
]

@dataclass
class ProteinRecord:
    uniprot_id: str
    gene_name: str
    protein_name: str
    subcellular_location: str
    sequence_length: Optional[float]
    idr_percentage: Optional[float]
    cc_percentage: Optional[float]
    row: Optional[Dict[str, Any]] = None


def get_db_connection():
    conn = getattr(g, "_db_conn", None)
    if conn is None:
        conn = g._db_conn = psycopg2.connect(SUPABASE_DB_URL)
    return conn


@app.teardown_appcontext
def close_db_connection(exception: Optional[BaseException]):
    conn = g.pop("_db_conn", None)
    if conn is not None:
        conn.close()


def normalise_row(row: Dict[str, Any], columns: Sequence[Tuple[str, str]]) -> Dict[str, Any]:
    normalized: Dict[str, Any] = {}
    for db_col, display_col in columns:
        normalized[display_col] = row.get(db_col)
    return normalized


def parse_page() -> int:
    try:
        value = int(request.args.get("page", 1))
    except (TypeError, ValueError):
        value = 1
    return max(1, value)


def get_page_size() -> int:
    try:
        per_page = int(request.args.get("per_page", DEFAULT_PAGE_SIZE))
    except (TypeError, ValueError):
        per_page = DEFAULT_PAGE_SIZE
    return max(1, min(per_page, MAX_PAGE_SIZE))


def parse_int_param(name: str, default: Optional[int] = None) -> Optional[int]:
    raw = request.args.get(name, "").strip()
    if raw == "":
        return default
    try:
        value = int(raw)
        return max(0, value)
    except ValueError:
        return default


def parse_float_param(name: str, default: Optional[float] = None) -> Optional[float]:
    raw = request.args.get(name, "").strip()
    if raw == "":
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def build_filter_conditions(
    search: Optional[str],
    search_mode: str,
    min_idr_length: Optional[int],
    min_cc_length: Optional[int],
    min_protein_length: Optional[int],
    max_protein_length: Optional[int],
    min_idr_pct: Optional[float],
    max_idr_pct: Optional[float],
    min_cc_pct: Optional[float],
    max_cc_pct: Optional[float],
    domain_term: Optional[str],
    location_term: Optional[str],
    hide_missing_protein: bool,
    alias: str = "p",
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []

    if search:
        term = f"%{search.lower()}%"
        mode = search_mode if search_mode in SEARCH_MODE_COLUMN_MAP else "all"
        columns = SEARCH_MODE_COLUMN_MAP.get(mode, SEARCH_MODE_COLUMN_MAP["all"])
        col_clauses = [f"LOWER({alias}.{col}) LIKE %s" for col in columns]
        clauses.append("(" + " OR ".join(col_clauses) + ")")
        params.extend([term] * len(col_clauses))

    if min_idr_length is not None:
        clauses.append(f"COALESCE({alias}.idr_residues, 0) >= %s")
        params.append(min_idr_length)

    if min_cc_length is not None:
        clauses.append(f"COALESCE({alias}.total_cc_length, 0) >= %s")
        params.append(min_cc_length)

    if min_protein_length is not None:
        clauses.append(f"COALESCE({alias}.sequence_length, 0) >= %s")
        params.append(min_protein_length)

    if max_protein_length is not None:
        clauses.append(f"COALESCE({alias}.sequence_length, 0) <= %s")
        params.append(max_protein_length)

    if min_idr_pct is not None:
        clauses.append(f"COALESCE({alias}.idr_percentage, 0) >= %s")
        params.append(min_idr_pct)
    if max_idr_pct is not None:
        clauses.append(f"COALESCE({alias}.idr_percentage, 0) <= %s")
        params.append(max_idr_pct)

    if min_cc_pct is not None:
        clauses.append(f"COALESCE({alias}.cc_percentage, 0) >= %s")
        params.append(min_cc_pct)
    if max_cc_pct is not None:
        clauses.append(f"COALESCE({alias}.cc_percentage, 0) <= %s")
        params.append(max_cc_pct)

    if domain_term:
        clauses.append(f"LOWER(COALESCE({alias}.domain_information, '')) LIKE %s")
        params.append(f"%{domain_term.lower()}%")

    if location_term:
        clauses.append(f"LOWER(COALESCE({alias}.subcellular_location, '')) LIKE %s")
        params.append(f"%{location_term.lower()}%")

    if hide_missing_protein:
        clauses.append(f"NULLIF(TRIM(COALESCE({alias}.protein_name, '')), '') IS NOT NULL")

    where_sql = "WHERE " + " AND ".join(clauses) if clauses else ""
    return where_sql, params


def build_idr_filter_conditions(
    search: Optional[str],
    min_len: Optional[int],
    max_len: Optional[int],
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []

    if search:
        term = f"%{search.lower()}%"
        clauses.append(
            "(" + " OR ".join(
                [
                    "LOWER(idr.uniprot_id) LIKE %s",
                    "LOWER(idr.gene_name) LIKE %s",
                    "LOWER(idr.protein_name) LIKE %s",
                ]
            ) + ")"
        )
        params.extend([term, term, term])

    if min_len is not None:
        clauses.append("COALESCE(idr.idr_length, 0) >= %s")
        params.append(min_len)
    if max_len is not None:
        clauses.append("COALESCE(idr.idr_length, 0) <= %s")
        params.append(max_len)

    where_sql = "WHERE " + " AND ".join(clauses) if clauses else ""
    return where_sql, params


def build_ppi_filter_conditions(
    source: Optional[str],
    uniprot_id: Optional[str],
    min_score: Optional[int],
    min_idr_pct: Optional[float],
    max_idr_pct: Optional[float],
    min_cc_pct: Optional[float],
    max_cc_pct: Optional[float],
    min_idr_len: Optional[int],
    min_cc_len: Optional[int],
    min_protein_len: Optional[int],
    max_protein_len: Optional[int],
    search: Optional[str],
    search_mode: str,
    hide_missing_protein: bool,
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []

    # Search across both partners
    if search:
        term = f"%{search.lower()}%"
        mode = search_mode if search_mode in SEARCH_MODE_COLUMN_MAP else "all"
        columns = SEARCH_MODE_COLUMN_MAP.get(mode, SEARCH_MODE_COLUMN_MAP["all"])
        col_clauses = []
        for col in columns:
            col_clauses.append(f"LOWER(p1.{col}) LIKE %s")
            col_clauses.append(f"LOWER(p2.{col}) LIKE %s")
        clauses.append("(" + " OR ".join(col_clauses) + ")")
        params.extend([term] * len(col_clauses))

    if source in {"biogrid", "string"}:
        clauses.append("e.source = %s")
        params.append(source)

    if uniprot_id:
        clauses.append("(e.uniprot_a = %s OR e.uniprot_b = %s)")
        params.extend([uniprot_id, uniprot_id])

    if min_score is not None:
        clauses.append("(e.source != 'string' OR COALESCE(e.combined_score, 0) >= %s)")
        params.append(min_score)

    # Protein length filters apply to both partners
    if min_protein_len is not None:
        clauses.append("COALESCE(p1.sequence_length, 0) >= %s")
        clauses.append("COALESCE(p2.sequence_length, 0) >= %s")
        params.extend([min_protein_len, min_protein_len])
    if max_protein_len is not None:
        clauses.append("COALESCE(p1.sequence_length, 0) <= %s")
        clauses.append("COALESCE(p2.sequence_length, 0) <= %s")
        params.extend([max_protein_len, max_protein_len])

    # Helper to build threshold parts
    def threshold_parts(alias: str, column: str, min_val: Optional[float], max_val: Optional[float]) -> Tuple[List[str], List[Any]]:
        parts: List[str] = []
        p: List[Any] = []
        if min_val is not None:
            parts.append(f"COALESCE({alias}.{column}, 0) >= %s")
            p.append(min_val)
        if max_val is not None:
            parts.append(f"COALESCE({alias}.{column}, 0) <= %s")
            p.append(max_val)
        return parts, p

    # IDR/CC percentage: one partner satisfies IDR range, the other CC range (either orientation)
    idr_pct_parts_p1, idr_pct_params_p1 = threshold_parts("p1", "idr_percentage", min_idr_pct, max_idr_pct)
    idr_pct_parts_p2, idr_pct_params_p2 = threshold_parts("p2", "idr_percentage", min_idr_pct, max_idr_pct)
    cc_pct_parts_p1, cc_pct_params_p1 = threshold_parts("p1", "cc_percentage", min_cc_pct, max_cc_pct)
    cc_pct_parts_p2, cc_pct_params_p2 = threshold_parts("p2", "cc_percentage", min_cc_pct, max_cc_pct)

    pct_orientation_clauses: List[str] = []
    pct_orientation_params: List[Any] = []
    orient1 = idr_pct_parts_p1 + cc_pct_parts_p2
    if orient1:
        pct_orientation_clauses.append("(" + " AND ".join(orient1) + ")")
        pct_orientation_params.extend(idr_pct_params_p1 + cc_pct_params_p2)
    orient2 = idr_pct_parts_p2 + cc_pct_parts_p1
    if orient2:
        pct_orientation_clauses.append("(" + " AND ".join(orient2) + ")")
        pct_orientation_params.extend(idr_pct_params_p2 + cc_pct_params_p1)
    if pct_orientation_clauses:
        clauses.append("(" + " OR ".join(pct_orientation_clauses) + ")")
        params.extend(pct_orientation_params)

    # IDR/CC length (aa): same orientation logic using idr_residues and total_cc_length
    idr_len_parts_p1, idr_len_params_p1 = threshold_parts("p1", "idr_residues", min_idr_len, None)
    idr_len_parts_p2, idr_len_params_p2 = threshold_parts("p2", "idr_residues", min_idr_len, None)
    cc_len_parts_p1, cc_len_params_p1 = threshold_parts("p1", "total_cc_length", min_cc_len, None)
    cc_len_parts_p2, cc_len_params_p2 = threshold_parts("p2", "total_cc_length", min_cc_len, None)

    len_orientation_clauses: List[str] = []
    len_orientation_params: List[Any] = []
    orient_len1 = idr_len_parts_p1 + cc_len_parts_p2
    if orient_len1:
        len_orientation_clauses.append("(" + " AND ".join(orient_len1) + ")")
        len_orientation_params.extend(idr_len_params_p1 + cc_len_params_p2)
    orient_len2 = idr_len_parts_p2 + cc_len_parts_p1
    if orient_len2:
        len_orientation_clauses.append("(" + " AND ".join(orient_len2) + ")")
        len_orientation_params.extend(idr_len_params_p2 + cc_len_params_p1)
    if len_orientation_clauses:
        clauses.append("(" + " OR ".join(len_orientation_clauses) + ")")
        params.extend(len_orientation_params)

    if hide_missing_protein:
        clauses.append("NULLIF(TRIM(COALESCE(p1.protein_name, '')), '') IS NOT NULL")
        clauses.append("NULLIF(TRIM(COALESCE(p2.protein_name, '')), '') IS NOT NULL")

    where_sql = "WHERE " + " AND ".join(clauses) if clauses else ""
    return where_sql, params


def paginate(total_items: int, page: int, per_page: int) -> Tuple[int, int, List[int], bool, bool, int]:
    total_pages = max(1, (total_items + per_page - 1) // per_page) if total_items else 1
    if total_items and page > total_pages:
        page = total_pages
    if total_items:
        start_item = (page - 1) * per_page + 1
        end_item = min(page * per_page, total_items)
        if total_pages <= 7:
            page_numbers = list(range(1, total_pages + 1))
            show_first = show_last = False
        else:
            start_page = max(1, page - 2)
            end_page = min(total_pages, page + 2)
            page_numbers = list(range(start_page, end_page + 1))
            show_first = start_page > 1
            show_last = end_page < total_pages
    else:
        start_item = end_item = 0
        page_numbers = []
        show_first = show_last = False
    return page, start_item, end_item, page_numbers, show_first, show_last, total_pages


def fetch_protein_page(
    table: str,
    filters: Tuple[Any, ...],
    page: int,
    per_page: int,
) -> Tuple[List[ProteinRecord], int, Tuple[str, List[Any]]]:
    (
        search,
        search_mode,
        min_idr_length,
        min_cc_length,
        min_protein_length,
        max_protein_length,
        min_idr_pct,
        max_idr_pct,
        min_cc_pct,
        max_cc_pct,
        domain_term,
        location_term,
        hide_missing,
    ) = filters
    where_sql, params = build_filter_conditions(
        search,
        search_mode,
        min_idr_length,
        min_cc_length,
        min_protein_length,
        max_protein_length,
        min_idr_pct,
        max_idr_pct,
        min_cc_pct,
        max_cc_pct,
        domain_term,
        location_term,
        hide_missing,
        alias="p",
    )
    conn = get_db_connection()
    total_items = 0
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM {table} p {where_sql}", params)
        total_items = cur.fetchone()[0]

    records: List[ProteinRecord] = []
    if total_items:
        offset = (page - 1) * per_page
        columns = [
            "uniprot_id",
            "gene_name",
            "protein_name",
            "subcellular_location",
            "sequence_length",
            "idr_percentage",
            "cc_percentage",
        ]
        select_sql = f"SELECT {', '.join(columns)} FROM {table} p {where_sql} ORDER BY p.ctid LIMIT %s OFFSET %s"
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(select_sql, params + [per_page, offset])
            for row in cur.fetchall():
                records.append(
                    ProteinRecord(
                        uniprot_id=row.get("uniprot_id", ""),
                        gene_name=row.get("gene_name", ""),
                        protein_name=row.get("protein_name", ""),
                        subcellular_location=row.get("subcellular_location", "") or "",
                        sequence_length=row.get("sequence_length"),
                        idr_percentage=row.get("idr_percentage"),
                        cc_percentage=row.get("cc_percentage"),
                    )
                )
    return records, total_items, (where_sql, params)


def _append_condition(where_clause: str, condition: str) -> str:
    if where_clause:
        return f"{where_clause} AND {condition}"
    return f"WHERE {condition}"


def compute_subcellular_counts(
    table: str,
    where_clause: str,
    params: Sequence[Any],
) -> Tuple[List[Tuple[str, int]], int, int]:
    conn = get_db_connection()
    base_params = tuple(params)
    total_sql = f"SELECT COUNT(*) FROM {table} p {where_clause}"
    unknown_condition = _append_condition(where_clause, "NULLIF(TRIM(COALESCE(p.subcellular_location, '')), '') IS NULL")
    counts_sql = f"""
        WITH filtered AS (
            SELECT p.uniprot_id, COALESCE(p.subcellular_location, '') AS subcellular_location
            FROM {table} p
            {where_clause}
        ),
        expanded AS (
            SELECT
                f.uniprot_id,
                NULLIF(TRIM(loc), '') AS location
            FROM filtered f
            CROSS JOIN LATERAL regexp_split_to_table(f.subcellular_location, ',') AS t(loc)
        ),
        dedup AS (
            SELECT DISTINCT
                uniprot_id,
                LOWER(location) AS normalized_location,
                location
            FROM expanded
            WHERE location IS NOT NULL
        )
        SELECT MIN(location) AS label, COUNT(*) AS freq
        FROM dedup
        GROUP BY normalized_location
        ORDER BY freq DESC
    """

    with conn.cursor() as cur:
        cur.execute(total_sql, base_params)
        total_entries = cur.fetchone()[0] or 0

        cur.execute(f"SELECT COUNT(*) FROM {table} p {unknown_condition}", base_params)
        unknown_count = cur.fetchone()[0] or 0

        cur.execute(counts_sql, base_params)
        rows = cur.fetchall()

    ordered = [(label, freq) for label, freq in rows if label]
    return ordered, total_entries, unknown_count


def fetch_protein_detail(uniprot_id: str) -> Optional[Dict[str, Any]]:
    conn = get_db_connection()
    for table in (PROTEINS_VER6, PROTEINS_VER9):
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                f"SELECT * FROM {table} WHERE UPPER(uniprot_id) = %s LIMIT 1",
                [uniprot_id.upper()],
            )
            row = cur.fetchone()
            if row:
                return normalise_row(row, PROTEIN_DISPLAY_COLUMNS)
    return None


def fetch_idr_detail(uniprot_id: str, number: int) -> Optional[Dict[str, Any]]:
    conn = get_db_connection()
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(
            f"SELECT * FROM {IDR_SEGMENTS_VER9} WHERE UPPER(uniprot_id) = %s AND idr_number = %s LIMIT 1",
            [uniprot_id.upper(), number],
        )
        row = cur.fetchone()
        if row:
            return normalise_row(row, IDR_DISPLAY_COLUMNS)
    return None


def extract_filters() -> Tuple[Any, ...]:
    search = request.args.get("search", "").strip()
    search_mode = request.args.get("search_mode", "all").strip().lower() or "all"
    if search_mode not in SEARCH_MODE_COLUMN_MAP:
        search_mode = "all"
    min_idr_length = parse_int_param("idr_min")
    min_cc_length = parse_int_param("cc_min")
    min_protein_length = parse_int_param("protein_len_min")
    max_protein_length = parse_int_param("protein_len_max")
    min_idr_pct = parse_float_param("idr_pct_min")
    max_idr_pct = parse_float_param("idr_pct_max")
    min_cc_pct = parse_float_param("cc_pct_min")
    max_cc_pct = parse_float_param("cc_pct_max")
    domain_term = request.args.get("domain_term", "").strip()
    location_term = request.args.get("location_term", "").strip()
    hide_missing_protein = request.args.get("hide_missing_protein", "").lower() in {"1", "true", "on"}
    return (
        search,
        search_mode,
        min_idr_length,
        min_cc_length,
        min_protein_length,
        max_protein_length,
        min_idr_pct,
        max_idr_pct,
        min_cc_pct,
        max_cc_pct,
        domain_term,
        location_term,
        hide_missing_protein,
    )


@app.context_processor
def inject_globals():
    return {
        "SEARCH_MODE_OPTIONS": SEARCH_MODE_OPTIONS,
        "ACCESS_PASSWORD_ENABLED": bool(ACCESS_PASSWORD),
    }


@app.before_request
def enforce_password_gate():
    if not ACCESS_PASSWORD:
        return
    if session.get("authenticated"):
        return
    endpoint = request.endpoint or ""
    if endpoint in {"login", "static", "health"}:
        return
    next_url = request.url
    return redirect(url_for("login", next=next_url))


def _safe_next_url(default: Optional[str] = None) -> str:
    candidate = request.args.get("next")
    if candidate:
        ref_url = urlparse(request.host_url)
        test_url = urlparse(urljoin(request.host_url, candidate))
        if test_url.scheme in {"http", "https"} and test_url.netloc == ref_url.netloc:
            return candidate
    return default or url_for("index")


def _safe_return_path(target: Optional[str], default: str) -> str:
    if not target:
        return default
    parsed = urlparse(target)
    if parsed.scheme or parsed.netloc:
        return default
    path = target if target.startswith("/") else f"/{target}"
    return path


@app.route("/login", methods=["GET", "POST"])
def login():
    error = None
    if request.method == "POST":
        password = request.form.get("password", "")
        if ACCESS_PASSWORD and password == ACCESS_PASSWORD:
            session["authenticated"] = True
            return redirect(_safe_next_url())
        error = "Invalid password."
    return render_template("login.html", error=error)


@app.route("/logout")
def logout():
    session.clear()
    return redirect(url_for("login"))


@app.route("/")
def index():
    page = parse_page()
    per_page = get_page_size()
    filters = extract_filters()
    (
        search,
        search_mode,
        idr_min,
        cc_min,
        protein_len_min,
        protein_len_max,
        idr_pct_min,
        idr_pct_max,
        cc_pct_min,
        cc_pct_max,
        domain_term,
        location_term,
        hide_missing_protein,
    ) = filters
    records, total_items, filter_clause = fetch_protein_page(PROTEINS_VER6, filters, page, per_page)
    page, start_item, end_item, page_numbers, show_first, show_last, total_pages = paginate(total_items, page, per_page)
    location_counts, location_total, unknown_count = compute_subcellular_counts(PROTEINS_VER6, *filter_clause)

    def page_url(target_page: int) -> str:
        return url_for(
            "index",
            page=target_page,
            per_page=per_page,
            search=search or None,
            search_mode=search_mode if search_mode != "all" else None,
            idr_min=idr_min if idr_min is not None else None,
            cc_min=cc_min if cc_min is not None else None,
            protein_len_min=protein_len_min if protein_len_min is not None else None,
            protein_len_max=protein_len_max if protein_len_max is not None else None,
            idr_pct_min=idr_pct_min if idr_pct_min is not None else None,
            idr_pct_max=idr_pct_max if idr_pct_max is not None else None,
            cc_pct_min=cc_pct_min if cc_pct_min is not None else None,
            cc_pct_max=cc_pct_max if cc_pct_max is not None else None,
            domain_term=domain_term or None,
            location_term=location_term or None,
            hide_missing_protein="1" if hide_missing_protein else None,
        )

    dataset_note = "Source: ver6 integrated dataset (Supabase)."
    current_path = request.full_path.rstrip("?") or request.path
    return render_template(
        "index.html",
        records=records,
        page=page,
        per_page=per_page,
        total_pages=total_pages,
        total_items=total_items,
        search=search,
        search_mode=search_mode,
        idr_min=idr_min,
        cc_min=cc_min,
        protein_len_min=protein_len_min,
        protein_len_max=protein_len_max,
        idr_pct_min=idr_pct_min,
        idr_pct_max=idr_pct_max,
        cc_pct_min=cc_pct_min,
        cc_pct_max=cc_pct_max,
        domain_term=domain_term,
        location_term=location_term,
        hide_missing_protein=hide_missing_protein,
        page_url=page_url,
        start_item=start_item,
        end_item=end_item,
        page_numbers=page_numbers,
        show_first=show_first,
        show_last=show_last,
        location_counts=location_counts,
        location_chart_url=url_for("subcellular_distribution_api"),
        dataset_note=dataset_note,
        page_title="Browse Entries",
        form_action=url_for("index"),
        location_chart_total=location_total,
        location_chart_unknown=unknown_count,
        current_path=current_path,
    )


@app.route("/canonical/")
def canonical_index():
    page = parse_page()
    per_page = get_page_size()
    filters = extract_filters()
    (
        search,
        search_mode,
        idr_min,
        cc_min,
        protein_len_min,
        protein_len_max,
        idr_pct_min,
        idr_pct_max,
        cc_pct_min,
        cc_pct_max,
        domain_term,
        location_term,
        hide_missing_protein,
    ) = filters
    records, total_items, filter_clause = fetch_protein_page(PROTEINS_VER9, filters, page, per_page)
    page, start_item, end_item, page_numbers, show_first, show_last, total_pages = paginate(total_items, page, per_page)
    location_counts, location_total, unknown_count = compute_subcellular_counts(PROTEINS_VER9, *filter_clause)

    def page_url(target_page: int) -> str:
        return url_for(
            "canonical_index",
            page=target_page,
            per_page=per_page,
            search=search or None,
            search_mode=search_mode if search_mode != "all" else None,
            idr_min=idr_min if idr_min is not None else None,
            cc_min=cc_min if cc_min is not None else None,
            protein_len_min=protein_len_min if protein_len_min is not None else None,
            protein_len_max=protein_len_max if protein_len_max is not None else None,
            idr_pct_min=idr_pct_min if idr_pct_min is not None else None,
            idr_pct_max=idr_pct_max if idr_pct_max is not None else None,
            cc_pct_min=cc_pct_min if cc_pct_min is not None else None,
            cc_pct_max=cc_pct_max if cc_pct_max is not None else None,
            domain_term=domain_term or None,
            location_term=location_term or None,
            hide_missing_protein="1" if hide_missing_protein else None,
        )

    dataset_note = "Canonical subset (ver9) · UniProt UP000005640."
    current_path = request.full_path.rstrip("?") or request.path
    return render_template(
        "index.html",
        records=records,
        page=page,
        per_page=per_page,
        total_pages=total_pages,
        total_items=total_items,
        search=search,
        search_mode=search_mode,
        idr_min=idr_min,
        cc_min=cc_min,
        protein_len_min=protein_len_min,
        protein_len_max=protein_len_max,
        idr_pct_min=idr_pct_min,
        idr_pct_max=idr_pct_max,
        cc_pct_min=cc_pct_min,
        cc_pct_max=cc_pct_max,
        domain_term=domain_term,
        location_term=location_term,
        hide_missing_protein=hide_missing_protein,
        page_url=page_url,
        start_item=start_item,
        end_item=end_item,
        page_numbers=page_numbers,
        show_first=show_first,
        show_last=show_last,
        location_counts=location_counts,
        location_chart_url=url_for("subcellular_distribution_api", dataset="canonical"),
        dataset_note=dataset_note,
        page_title="Canonical Entries",
        form_action=url_for("canonical_index"),
        location_chart_total=location_total,
        location_chart_unknown=unknown_count,
        current_path=current_path,
    )


def get_threshold_plot(uniprot_id: str) -> Optional[str]:
    candidates = {uniprot_id.upper()}
    raw = uniprot_id
    for delim in [";", ",", " "]:
        if delim in raw:
            parts = [part.strip().upper() for part in raw.split(delim)]
            candidates.update(part for part in parts if part)
    for candidate in candidates:
        plot = threshold_images.get(candidate)
        if plot:
            return plot
    return None


@app.route("/protein/<uniprot_id>")
def protein_detail(uniprot_id: str):
    row = fetch_protein_detail(uniprot_id)
    if not row:
        abort(404)
    record = ProteinRecord(
        uniprot_id=row.get("UniProt_ID", ""),
        gene_name=row.get("Gene_Name", ""),
        protein_name=row.get("Protein_Name", ""),
        subcellular_location=row.get("Subcellular_Location", "") or "",
        sequence_length=row.get("Sequence_Length"),
        idr_percentage=row.get("IDR_Percentage"),
        cc_percentage=row.get("CC_Percentage"),
        row=row,
    )
    plot_filename = get_threshold_plot(uniprot_id)
    plot_url = url_for("analysis_plot", filename=plot_filename) if plot_filename else None
    return_to = _safe_return_path(request.args.get("return_to"), url_for("index"))
    return render_template("detail.html", record=record, plot_url=plot_url, return_to=return_to)


@app.route("/analysis/<path:filename>")
def analysis_plot(filename: str):
    if not ANALYSIS_DIR.exists():
        abort(404)
    return send_from_directory(ANALYSIS_DIR, filename, as_attachment=False)


@app.route("/health")
def health():
    return ("ok", 200, {"Content-Type": "text/plain; charset=utf-8"})


@app.route("/api/subcellular_distribution")
def subcellular_distribution_api():
    dataset = request.args.get("dataset", "").strip().lower()
    table = PROTEINS_VER9 if dataset == "canonical" else PROTEINS_VER6
    filters = extract_filters()
    where_sql, params = build_filter_conditions(*filters, alias="p")
    counts, total_entries, unknown_count = compute_subcellular_counts(table, where_sql, params)
    labels = [label for label, _ in counts]
    values = [value for _, value in counts]
    known_total = max(total_entries - unknown_count, 0)
    return jsonify(
        {
            "labels": labels,
            "values": values,
            "total": known_total,
            "total_known": known_total,
            "total_entries": total_entries,
            "unknown_count": unknown_count,
        }
    )


@app.route("/supramolecular")
def supramolecular():
    page = parse_page()
    per_page = get_page_size()
    search = request.args.get("search", "").strip()
    search_mode = request.args.get("search_mode", "all").strip().lower() or "all"
    if search_mode not in SEARCH_MODE_COLUMN_MAP:
        search_mode = "all"
    source = request.args.get("source", "").strip().lower()
    source = source if source in {"biogrid", "string"} else None
    uniprot_id = request.args.get("uniprot", "").strip().upper() or None
    min_score = parse_int_param("score_min")
    min_idr_len = parse_int_param("idr_len_min")
    min_cc_len = parse_int_param("cc_len_min")
    min_idr_pct = parse_float_param("idr_pct_min")
    max_idr_pct = parse_float_param("idr_pct_max")
    min_cc_pct = parse_float_param("cc_pct_min")
    max_cc_pct = parse_float_param("cc_pct_max")
    min_protein_len = parse_int_param("protein_len_min")
    max_protein_len = parse_int_param("protein_len_max")
    hide_missing_protein = request.args.get("hide_missing_protein", "").lower() in {"1", "true", "on"}

    where_sql, params = build_ppi_filter_conditions(
        source,
        uniprot_id,
        min_score,
        min_idr_pct,
        max_idr_pct,
        min_cc_pct,
        max_cc_pct,
        min_idr_len,
        min_cc_len,
        min_protein_len,
        max_protein_len,
        search,
        search_mode,
        hide_missing_protein,
    )

    conn = get_db_connection()
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT COUNT(*)
            FROM {PPI_EDGES} e
            JOIN {PROTEINS_VER10} p1 ON e.uniprot_a = p1.uniprot_id
            JOIN {PROTEINS_VER10} p2 ON e.uniprot_b = p2.uniprot_id
            {where_sql}
            """,
            params,
        )
        total_items = cur.fetchone()[0]

    records: List[Dict[str, Any]] = []
    if total_items:
        offset = (page - 1) * per_page
        query = f"""
            SELECT
                e.source,
                e.combined_score,
                p1.uniprot_id AS a_id,
                p1.gene_name AS a_gene,
                p1.idr_percentage AS a_idr_pct,
                p1.cc_percentage AS a_cc_pct,
                p2.uniprot_id AS b_id,
                p2.gene_name AS b_gene,
                p2.idr_percentage AS b_idr_pct,
                p2.cc_percentage AS b_cc_pct
            FROM {PPI_EDGES} e
            JOIN {PROTEINS_VER10} p1 ON e.uniprot_a = p1.uniprot_id
            JOIN {PROTEINS_VER10} p2 ON e.uniprot_b = p2.uniprot_id
            {where_sql}
            ORDER BY e.source, e.combined_score DESC NULLS LAST, p1.uniprot_id, p2.uniprot_id
            LIMIT %s OFFSET %s
        """
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(query, params + [per_page, offset])
            rows = cur.fetchall()
            for r in rows:
                records.append(
                    {
                        "source": r.get("source"),
                        "combined_score": r.get("combined_score"),
                        "a_id": r.get("a_id"),
                        "a_gene": r.get("a_gene"),
                        "a_idr_pct": r.get("a_idr_pct"),
                        "a_cc_pct": r.get("a_cc_pct"),
                        "b_id": r.get("b_id"),
                        "b_gene": r.get("b_gene"),
                        "b_idr_pct": r.get("b_idr_pct"),
                        "b_cc_pct": r.get("b_cc_pct"),
                    }
                )

    page, start_item, end_item, page_numbers, show_first, show_last, total_pages = paginate(total_items, page, per_page)

    def page_url(target_page: int) -> str:
        return url_for(
            "supramolecular",
            page=target_page,
            per_page=per_page,
            source=source or None,
            uniprot=uniprot_id or None,
            score_min=min_score if min_score is not None else None,
            idr_pct_min=min_idr_pct if min_idr_pct is not None else None,
            cc_pct_min=min_cc_pct if min_cc_pct is not None else None,
        )

    return render_template(
        "supramolecular.html",
        records=records,
        total_items=total_items,
        page=page,
        start_item=start_item,
        end_item=end_item,
        page_numbers=page_numbers,
        show_first=show_first,
        show_last=show_last,
        total_pages=total_pages,
        per_page=per_page,
        source=source,
        uniprot_id=uniprot_id,
        min_score=min_score,
        min_idr_len=min_idr_len,
        min_cc_len=min_cc_len,
        min_idr_pct=min_idr_pct,
        max_idr_pct=max_idr_pct,
        min_cc_pct=min_cc_pct,
        max_cc_pct=max_cc_pct,
        min_protein_len=min_protein_len,
        max_protein_len=max_protein_len,
        search=search,
        search_mode=search_mode,
        hide_missing_protein=hide_missing_protein,
        page_url=page_url,
    )


@app.route("/idr/")
def idr_index():
    page = parse_page()
    per_page = get_page_size()
    search = request.args.get("search", "").strip()
    min_len = parse_int_param("length_min")
    max_len = parse_int_param("length_max")
    where_sql, params = build_idr_filter_conditions(search, min_len, max_len)

    conn = get_db_connection()
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM {IDR_SEGMENTS_VER9} idr {where_sql}", params)
        total_items = cur.fetchone()[0]

    records: List[Dict[str, Any]] = []
    if total_items:
        offset = (page - 1) * per_page
        query = f"""
            SELECT *
            FROM {IDR_SEGMENTS_VER9} idr
            {where_sql}
            ORDER BY idr.ctid
            LIMIT %s OFFSET %s
        """
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(query, params + [per_page, offset])
            rows = cur.fetchall()
            records = [normalise_row(row, IDR_DISPLAY_COLUMNS) for row in rows]

    page, start_item, end_item, page_numbers, show_first, show_last, total_pages = paginate(total_items, page, per_page)

    def page_url(target_page: int) -> str:
        return url_for(
            "idr_index",
            page=target_page,
            per_page=per_page,
            search=search or None,
            length_min=min_len if min_len is not None else None,
            length_max=max_len if max_len is not None else None,
        )

    dataset_note = "IDR segments (≥30 aa) from canonical ver9."
    current_path = request.full_path.rstrip("?") or request.path
    return render_template(
        "idr_index.html",
        df=records,
        page=page,
        per_page=per_page,
        total_pages=total_pages,
        total_items=total_items,
        start_item=start_item,
        end_item=end_item,
        page_url=page_url,
        page_numbers=page_numbers,
        show_first=show_first,
        show_last=show_last,
        search=search,
        length_min=min_len,
        length_max=max_len,
        dataset_note=dataset_note,
        current_path=current_path,
    )


@app.route("/idr/<uid>/<int:number>")
def idr_detail(uid: str, number: int):
    row = fetch_idr_detail(uid, number)
    if not row:
        abort(404)
    return render_template(
        "idr_detail.html",
        record=row,
        tsne_image_available=False,
        tsne_marker=None,
        return_to=_safe_return_path(request.args.get("return_to"), url_for("idr_index")),
    )


threshold_images: Dict[str, str] = {}
if ANALYSIS_DIR.exists():
    ids = sorted(
        {
            path.name.split("_")[0].upper()
            for path in ANALYSIS_DIR.glob("*_combined_analysis_threshold.png")
        }
    )
    for uid in ids[:100]:
        threshold_images[uid] = f"{uid}_combined_analysis_threshold.png"


if __name__ == "__main__":
    app.run(debug=True)
