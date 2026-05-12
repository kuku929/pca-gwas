"""
msprime PCA failure simulations.

Run all experiments:   uv run main.py
Run one experiment:    uv run main.py timing | ghost | stepping | flow
"""
import sys
from pathlib import Path

IDEAS = {
    "timing":    "ideas.admixture_timing",
    "ghost":     "ideas.ghost_population",
    "ghost_sweep": "ideas.ghost_c_sweep",
    "stepping":  "ideas.stepping_stone",
    "flow":      "ideas.continuous_flow",
}


def main():
    targets = sys.argv[1:] or list(IDEAS)
    unknown = [t for t in targets if t not in IDEAS]
    if unknown:
        print(f"Unknown experiments: {unknown}")
        print(f"Valid: {list(IDEAS)}")
        sys.exit(1)

    results = {}
    for name in targets:
        import importlib
        mod = importlib.import_module(IDEAS[name])
        results[name] = mod.run()

    print("\n========== All done ==========")
    for name, lams in results.items():
        print(f"\n[{name}]")
        if isinstance(lams, dict):
            for scenario, vals in lams.items():
                if isinstance(vals, dict):
                    print(f"  {scenario}: {vals}")
                else:
                    print(f"  {vals}")


if __name__ == "__main__":
    main()
