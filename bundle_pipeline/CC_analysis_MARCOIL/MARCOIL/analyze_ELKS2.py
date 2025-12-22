import re
import matplotlib.pyplot as plt
import pandas as pd

# Read ProbList file
with open('MARCOIL/Outputs/ProbList', 'r') as f:
    lines = f.readlines()

# Extract probability data
data = []
for line in lines:
    match = re.match(r'\s*(\d+)\s+(\w)\s+([\d.]+)\s+(\w)', line)
    if match:
        position = int(match.group(1))
        residue = match.group(2)
        probability = float(match.group(3))
        heptad_phase = match.group(4)
        data.append({
            'Position': position,
            'Residue': residue,
            'Probability': probability,
            'Heptad_Phase': heptad_phase
        })

# Create DataFrame
df = pd.DataFrame(data)

# Save to CSV
df.to_csv('ELKS2_SS306_probabilities.csv', index=False)
print(f"Saved numerical data to ELKS2_SS306_probabilities.csv ({len(df)} positions)")

# Create visualization
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

# Plot 1: Coiled-coil probability
ax1.plot(df['Position'], df['Probability'], linewidth=1.5, color='blue')
ax1.fill_between(df['Position'], df['Probability'], alpha=0.3, color='blue')
ax1.axhline(y=50, color='red', linestyle='--', linewidth=1, label='50% threshold')
ax1.set_xlabel('Residue Position', fontsize=12)
ax1.set_ylabel('Coiled-Coil Probability (%)', fontsize=12)
ax1.set_title('ELKS2 (SS306) Coiled-Coil Prediction (MARCOIL)', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_ylim(0, 105)

# Plot 2: High confidence regions (>80%)
high_conf = df[df['Probability'] >= 80].copy()
ax2.scatter(high_conf['Position'], high_conf['Probability'], 
           c=high_conf['Probability'], cmap='YlOrRd', s=20, alpha=0.6)
ax2.axhline(y=80, color='orange', linestyle='--', linewidth=1, label='80% threshold')
ax2.axhline(y=100, color='red', linestyle='--', linewidth=1, label='100%')
ax2.set_xlabel('Residue Position', fontsize=12)
ax2.set_ylabel('Coiled-Coil Probability (%)', fontsize=12)
ax2.set_title('High Confidence Coiled-Coil Regions (≥80%)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_ylim(75, 105)

plt.tight_layout()
plt.savefig('ELKS2_SS306_coiled_coil_analysis.png', dpi=300, bbox_inches='tight')
print("Saved graph to ELKS2_SS306_coiled_coil_analysis.png")

# Generate summary statistics
print("\n=== Summary Statistics ===")
print(f"Total residues analyzed: {len(df)}")
print(f"Residues with >50% probability: {len(df[df['Probability'] > 50])} ({len(df[df['Probability'] > 50])/len(df)*100:.1f}%)")
print(f"Residues with >80% probability: {len(df[df['Probability'] > 80])} ({len(df[df['Probability'] > 80])/len(df)*100:.1f}%)")
print(f"Residues with 100% probability: {len(df[df['Probability'] == 100])} ({len(df[df['Probability'] == 100])/len(df)*100:.1f}%)")
print(f"\nAverage probability: {df['Probability'].mean():.2f}%")
print(f"Maximum probability: {df['Probability'].max():.2f}%")
print(f"Minimum probability: {df['Probability'].min():.2f}%")

# Identify continuous coiled-coil regions
regions = []
in_region = False
start = None
for i, row in df.iterrows():
    if row['Probability'] >= 80:
        if not in_region:
            start = row['Position']
            in_region = True
    else:
        if in_region:
            regions.append((start, df.iloc[i-1]['Position']))
            in_region = False
if in_region:
    regions.append((start, df.iloc[-1]['Position']))

print(f"\n=== Predicted Coiled-Coil Domains (≥80% probability) ===")
for i, (s, e) in enumerate(regions, 1):
    print(f"Region {i}: residues {s}-{e} (length: {e-s+1})")
