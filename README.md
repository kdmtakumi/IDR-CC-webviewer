# IDR+CC Database Web Viewer

This is a lightweight Flask application for browsing the merged IDR+CC dataset stored on Supabase. It exposes the ver6 integrated table (`proteins_ver6`), the canonical ver9 table (`proteins_ver9`), and the canonical IDR segment list (`idr_segments_ver9`). The UI supports keyword filtering, pagination, detailed drill-down pages, an IDR-specific browser, and subcellular-location charts.

## Structure
- `app.py`: Flask entry point and routes
- `templates/`: DisProt-inspired UI templates (`layout.html`, `index.html`, `detail.html`, `idr_index.html`, `idr_detail.html`, `tsne.html`)
- `static/style.css`: Styling for all pages
- `requirements.txt`: Minimal dependencies (Flask, psycopg2-binary, chart.js via CDN)
- `.venv/`: Optional local virtual environment

## Getting Started
```bash
cd "/Users/kodamatakumi/Desktop/IDR+CC DB/webviewer"
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
export SUPABASE_DB_URL="postgresql://<user>:<password>@<host>:<port>/<db>?sslmode=require"
export IDRCC_PASSWORD="ShimoLAB0501"        # optional gate
FLASK_APP=app.py flask run
```
Open `http://127.0.0.1:5000/` and authenticate with the shared password (default `ShimoLAB0501`).

## Features
- Keyword search (UniProt ID / Gene Name / Protein Name / Subcellular Location) with selectable search scope
- Pagination controls (10, 25, 50, 100 rows)
- Filters for IDR/CC absolute length and percentage, domain/location keywords, and “hide missing protein names”
- Protein detail pages with threshold plots (when available)
- Dedicated IDR list (≥30 aa) with cluster/distance metadata and per-IDR detail view
- Subcellular location pie chart that respects current filters
- Canonical view (ver9) covering all 83,607 canonical entries

## Deploying on Render + Supabase
1. **Supabase**
   - Create the `proteins_ver6`, `proteins_ver9`, and `idr_segments_ver9` tables (see `app.py` for column lists).
   - Upload the CSV files via `\copy` or Supabase Studio.
2. **Render**
   - Connect the GitHub repository and create a “Web Service”.
   - Build command: `pip install --upgrade pip && pip install -r requirements.txt`
   - Start command: `gunicorn app:app --chdir webviewer --bind 0.0.0.0:$PORT`
   - Environment variables:
     - `SUPABASE_DB_URL` (connection string, preferably using the pooling host)
     - `IDRCC_PASSWORD` (shared password) and `IDRCC_SECRET_KEY`
3. Set the Render instance to the same region as Supabase (or nearby) for lower latency.

### Health Check (Render / UptimeRobot)
- Purpose: external uptime monitoring without authentication; `/health` returns plain `ok` and is exempt from the login gate.
- UptimeRobot example: HTTP(s) monitor → URL `https://idr-cc-webviewer.onrender.com/health` → interval 5 minutes → Custom HTTP Header `User-Agent: Mozilla/5.0` (任意の一般的な UA) → 通知先のみ選択。
- Access log負荷を抑えるため、監視間隔は5分程度を推奨（現状のログ量は微小）。
- Renderログで確認する場合: ダッシュボードの Logs を開き、`/health` でフィルタ。5分刻みで 200 応答が並んでいれば正常。

## Data Notes
- Both “Browse” (ver6) and Canonical (ver9) tables originate from UniProt proteome UP000005640.
- The canonical FASTA (`uniprotkb_proteome_UP000005640_2025_10_29_canonical.fasta`) plus ver6 annotations were merged to create `all_human_protein_database_with_IDR-CCinformation_ver9.csv` (83,607 entries).
- The additional 3,070 entries introduced in ver9 are stored in `20251029_additional_analysis_filled.csv` (source snapshots in `20251029_additional_analysis.csv`).
- `IDR_sequences_30aa+_ver9.csv` interprets canonical `IDR_Boundaries` as 0-based half-open intervals, extracts segments ≥30 aa, preserves all metadata, and assigns sequential `IDR_Number` values within each protein.

## Outstanding Tasks
- Some of the 3,070 additional entries still lack gene/protein/location/domain annotations because the information is not yet released on UniProt. Fetch the latest metadata through the UniProt API when needed.

## Deployment (Render)
1. Push your code to GitHub (`kdmtakumi/IDR-CC-webviewer`).
2. On Render, create a new **Web Service** with Runtime = Native.
3. Build command: `pip install --upgrade pip && pip install -r requirements.txt`
4. Start command: `gunicorn app:app --chdir webviewer --bind 0.0.0.0:$PORT`
5. Environment variables:
   - `SUPABASE_DB_URL` = `postgresql://postgres.stjtqqlqpoxcpzywrxrz:ShimoLAB0501@aws-1-ap-southeast-1.pooler.supabase.com:6543/postgres?sslmode=require&client_encoding=utf8`
   - `IDRCC_PASSWORD` = `ShimoLAB0501`
   - `IDRCC_SECRET_KEY` = (任意の長い文字列)
6. Deploy and access https://idr-cc-webviewer.onrender.com
