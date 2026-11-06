"""Behaviour scenarios: small, simulated plasmids that each demonstrate one choice or trade-off.

Each scenario builds a plasmid and its contigs from random sequence (seeded), runs the real
``plasmid-surgery`` command, and returns what it found. ``tests/test_scenarios.py`` checks every claim;
``docs/gen_report.py`` turns the same runs into the documentation's Behaviour report, with figures. So
the report can't drift from what the code does: if a behaviour changes, a test fails.
"""

from tests.scenarios import assembly, features, find_refs, indels, origin, overlaps, reference_choice, repeats

MODULES = (origin, repeats, overlaps, indels, assembly, features, reference_choice, find_refs)

# In report order.
SCENARIOS = [m.SCENARIO for m in MODULES]
