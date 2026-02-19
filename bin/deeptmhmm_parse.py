import csv
import re
import sys

input_file = sys.argv[1]  # e.g., 'biolib_results/predicted_topologies.3line'
output_file = sys.argv[2]  # e.g., 'deeptmhmm_summary.txt'

def find_tm_helices(topology):
    """Find TM helices as stretches of 'M'"""
    helices = []
    for match in re.finditer(r"M+", topology):
        start = match.start() + 1
        end = match.end()
        helices.append((start, end))
    return helices

with open(input_file) as f, open(output_file, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "Sequence_ID",
        "Length",
        "Num_TM_helices",
        "TM_helices(start-end)",
        "Prediction",
        "Topology"
    ])

    lines = f.readlines()
    for i in range(0, len(lines), 3):
        # The first line contains ID and prediction label separated by ' | '
        header = lines[i].strip()
        if " | " in header:
            seq_id, prediction = header[1:].split(" | ")
        else:
            seq_id = header[1:]
            prediction = ""
        seq = lines[i+1].strip()
        topology = lines[i+2].strip()

        tm_coords = find_tm_helices(topology)
        tm_count = len(tm_coords)
        tm_str = ";".join([f"{s}-{e}" for s, e in tm_coords])

        writer.writerow([
            seq_id,
            len(seq),
            tm_count,
            tm_str,
            prediction,
            topology
        ])

print(f"Saved summary to {output_file}")
