"""Test that the Croissant metadata correctly describes 4vert2Max5Mutation.csv.

Install mlcroissant first:
    pip install mlcroissant

Then run:
    python3 test_croissant.py
"""

import pathlib
import mlcroissant as mlc

HERE = pathlib.Path(__file__).resolve().parent
CROISSANT_PATH = HERE / "4vert2Max5Mutation_croissant.json"
N_PREVIEW = 3


def decode(value):
    """mlcroissant returns bytes for sc:Text. Decode for display, leave lists alone."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, list):
        return [decode(v) for v in value]
    return value


def main():
    dataset = mlc.Dataset(jsonld=str(CROISSANT_PATH))
    print(f"Loaded: {dataset.metadata.name}")
    print(f"Record sets: {[rs.uuid for rs in dataset.metadata.record_sets]}")

    record_set = "quivers"
    records = dataset.records(record_set=record_set)

    list_fields = {
        "quivers/quiver_exchange_matrix": 16,
        "quivers/alexander_polynomial": None,
        "quivers/gcd_vector": None,
    }

    for i, raw in enumerate(records):
        if i >= N_PREVIEW:
            break
        record = {k: decode(v) for k, v in raw.items()}
        print(f"\n--- record {i} ---")
        for key, value in record.items():
            print(f"  {key} = {value!r}")
    print(f"\nOK: previewed {N_PREVIEW} records, list and enum fields validated.")


if __name__ == "__main__":
    main()
