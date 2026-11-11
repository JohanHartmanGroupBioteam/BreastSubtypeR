"""Independent decimal-token audit; writes only diagnostic artifacts."""
import csv
import decimal
import json
import math
import pathlib
import struct
import sys

source, destination = map(pathlib.Path, sys.argv[1:3])
decimal.getcontext().prec = 1200
with source.open(newline="", encoding="utf-8-sig") as stream:
    table = list(csv.reader(stream, delimiter="\t"))
header, records = table[0], table[1:]
nrow, ncol = len(records), len(records[0]) - 1
assert all(len(row) == ncol + 1 for row in records)
assert len(header) == ncol + 1
arrays = {}
for name in ("reference", "read", "parsed"):
    data = (destination / (name + ".bin")).read_bytes()
    assert len(data) == 8 * nrow * ncol
    arrays[name] = struct.unpack("<" + "d" * (nrow * ncol), data)
counts = dict(cells=nrow * ncol, nearest_proof_failed=0, reference_differs=0,
              read_differs=0, parser_differs=0, read_parser_differs=0)
nearest_bytes = bytearray()
max_reference_decimal_error = decimal.Decimal(0)
max_read_decimal_error = decimal.Decimal(0)
with (destination / "decimal-differences.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(("index_R", "row", "column", "probe", "sample", "token",
                     "nearest_hex", "reference_hex", "read_hex", "parser_hex",
                     "token_minus_reference", "token_minus_read", "nearest_proved"))
    for row, record in enumerate(records):
        for column, token in enumerate(record[1:]):
            exact = decimal.Decimal(token)
            nearest = float(token)
            assert math.isfinite(nearest)
            nearest_bytes.extend(struct.pack("<d", nearest))
            index = column * nrow + row  # R's column-major storage
            reference, read, parsed = (arrays[k][index] for k in ("reference", "read", "parsed"))
            exact_nearest = decimal.Decimal.from_float(nearest)
            previous = decimal.Decimal.from_float(math.nextafter(nearest, -math.inf))
            following = decimal.Decimal.from_float(math.nextafter(nearest, math.inf))
            distance = abs(exact - exact_nearest)
            neighbors = (abs(exact - previous), abs(exact - following))
            even = struct.unpack("<Q", struct.pack("<d", nearest))[0] % 2 == 0
            proved = all(distance < d or (distance == d and even) for d in neighbors)
            counts["nearest_proof_failed"] += not proved
            counts["reference_differs"] += nearest != reference
            counts["read_differs"] += nearest != read
            counts["parser_differs"] += nearest != parsed
            counts["read_parser_differs"] += read != parsed
            ref_error = exact - decimal.Decimal.from_float(reference)
            read_error = exact - decimal.Decimal.from_float(read)
            max_reference_decimal_error = max(max_reference_decimal_error, abs(ref_error))
            max_read_decimal_error = max(max_read_decimal_error, abs(read_error))
            if nearest != reference or nearest != read or nearest != parsed or not proved:
                writer.writerow((index + 1, row + 1, column + 1, record[0], header[column + 1],
                                 token, nearest.hex(), reference.hex(), read.hex(), parsed.hex(),
                                 str(ref_error), str(read_error), proved))
(destination / "nearest.bin").write_bytes(nearest_bytes)
counts.update(python=sys.version, rows=nrow, columns=ncol, decimal_precision=1200,
              max_reference_decimal_error=str(max_reference_decimal_error),
              max_read_decimal_error=str(max_read_decimal_error))
(destination / "decimal-summary.json").write_text(json.dumps(counts, indent=2) + "\n", encoding="utf-8")
assert counts["nearest_proof_failed"] == 0
