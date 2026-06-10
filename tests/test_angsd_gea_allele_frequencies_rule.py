from pathlib import Path


RULE_FILE = Path(__file__).resolve().parents[1] / "workflow" / "rules" / "angsd_gea_allele_frequencies.smk"


def test_all_samples_angsd_command_uses_large_dataset_safe_flags():
    rule_text = RULE_FILE.read_text()
    first_rule = rule_text.split("rule create_sites_file:", 1)[0]

    assert "-nQueueSize 50" in first_rule
    assert "-dosnpstat 1" not in first_rule
