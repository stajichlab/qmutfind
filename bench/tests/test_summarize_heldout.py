import math

from summarize_heldout import summarize, wilson


def test_wilson_interval_matches_a_hand_computed_value_and_handles_the_edges():
    lo, hi = wilson(8, 10)
    assert abs(lo - 0.4901) < 1e-3 and abs(hi - 0.9433) < 1e-3

    lo, hi = wilson(0, 10)
    assert lo == 0.0 and abs(hi - 0.2775) < 1e-3

    lo, hi = wilson(10, 10)
    assert abs(lo - 0.7225) < 1e-3 and hi == 1.0

    assert all(math.isnan(x) for x in wilson(0, 0))


ROWS = [
    # seed, mode, tech, coverage, klass, truth, tp, calls, fp
    dict(seed=3, mode="standard", tech="art", coverage=30, klass="SNP", truth=100, tp=90, calls=92, fp=2),
    dict(seed=4, mode="standard", tech="art", coverage=30, klass="SNP", truth=100, tp=100, calls=100, fp=0),
    dict(seed=3, mode="standard", tech="art", coverage=10, klass="SNP", truth=100, tp=50, calls=50, fp=0),
]


def test_results_are_pooled_over_seeds_with_intervals_and_the_per_seed_range():
    table = summarize(ROWS)

    row = next(r for r in table if r["coverage"] == 30)
    assert (row["truth"], row["tp"], row["calls"], row["fp"], row["seeds"]) == (200, 190, 192, 2, 2)
    assert abs(row["sensitivity"] - 0.95) < 1e-9
    assert abs(row["precision"] - 190 / 192) < 1e-9
    assert (row["sens_min"], row["sens_max"]) == (0.9, 1.0)
    assert row["sens_lo"] < 0.95 < row["sens_hi"]


def test_rows_are_grouped_by_mode_technology_coverage_and_class_and_sorted():
    table = summarize(ROWS)

    assert [(r["mode"], r["tech"], r["coverage"], r["klass"]) for r in table] == [
        ("standard", "art", 10, "SNP"),
        ("standard", "art", 30, "SNP"),
    ]


def test_precision_is_nan_when_there_are_no_calls():
    table = summarize([dict(seed=3, mode="standard", tech="ont", coverage=6, klass="LARGE_INDEL", truth=10, tp=0, calls=0, fp=0)])

    assert math.isnan(table[0]["precision"]) and table[0]["sensitivity"] == 0.0
