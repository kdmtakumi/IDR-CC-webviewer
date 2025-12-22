import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Read ProbList file
with open('MARCOIL/Outputs/ProbList', 'r') as f:
    lines = f.readlines()

# Extract probability data
positions = []
residues = []
probabilities = []
heptad_phases = []

for line in lines:
    match = re.match(r'\s*(\d+)\s+(\w)\s+([\d.]+)\s+(\w)', line)
    if match:
        positions.append(int(match.group(1)))
        residues.append(match.group(2))
        probabilities.append(float(match.group(3)))
        heptad_phases.append(match.group(4))

# Save to CSV
with open('ELKS2_SS306_probabilities.csv', 'w') as f:
    f.write('Position,Residue,Probability,Heptad_Phase\n')
    for i in range(len(positions)):
        f.write(f'{positions[i]},{residues[i]},{probabilities[i]},{heptad_phases[i]}\n')

print(f"Saved numerical data to ELKS2_SS306_probabilities.csv ({len(positions)} positions)")

# Create visualization
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

# Plot 1: Coiled-coil probability
ax1.plot(positions, probabilities, linewidth=1.5, color='blue')
ax1.fill_between(positions, probabilities, alpha=0.3, color='blue')
ax1.axhline(y=50, color='red', linestyle='--', linewidth=1, label='50% threshold')
ax1.set_xlabel('Residue Position', fontsize=12)
ax1.set_ylabel('Coiled-Coil Probability (%)', fontsize=12)
ax1.set_title('ELKS2 (SS306) Coiled-Coil Prediction (MARCOIL)', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_ylim(0, 105)

# Plot 2: High confidence regions (>80%)
high_conf_pos = [positions[i] for i in range(len(positions)) if probabilities[i] >= 80]
high_conf_prob = [probabilities[i] for i in range(len(positions)) if probabilities[i] >= 80]
ax2.scatter(high_conf_pos, high_conf_prob, c=high_conf_prob, cmap='YlOrRd', s=20, alpha=0.6)
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
total = len(probabilities)
over50 = len([p for p in probabilities if p > 50])
over80 = len([p for p in probabilities if p > 80])
perfect = len([p for p in probabilities if p == 100])

print("\n=== Summary Statistics ===")
print(f"Total residues analyzed: {total}")
print(f"Residues with >50% probability: {over50} ({over50/total*100:.1f}%)")
print(f"Residues with >80% probability: {over80} ({over80/total*100:.1f}%)")
print(f"Residues with 100% probability: {perfect} ({perfect/total*100:.1f}%)")
print(f"\nAverage probability: {sum(probabilities)/len(probabilities):.2f}%")
print(f"Maximum probability: {max(probabilities):.2f}%")
print(f"Minimum probability: {min(probabilities):.2f}%")

# Identify continuous coiled-coil regions
regions = []
in_region = False
start = None
for i in range(len(positions)):
    if probabilities[i] >= 80:
        if not in_region:
            start = positions[i]
            in_region = True
    else:
        if in_region:
            regions.append((start, positions[i-1]))
            in_region = False
if in_region:
    regions.append((start, positions[-1]))

print(f"\n=== Predicted Coiled-Coil Domains (≥80% probability) ===")
for i, (s, e) in enumerate(regions, 1):
    print(f"Region {i}: residues {s}-{e} (length: {e-s+1})")
