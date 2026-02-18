#!/usr/bin/env python3
"""Test script to verify the metalotl app loads correctly."""

try:
    print("Testing metalotl import...")
    from metalotl.app import app
    print(f"✓ App loaded successfully! Type: {type(app)}")
    
    print("\nTesting constants...")
    from metalotl._constants import DATA, GENES
    print(f"✓ Samples loaded: {len(DATA)}")
    print(f"✓ Genes loaded: {len(GENES)}")
    print(f"\nSamples:")
    for s in DATA:
        print(f"  - {s}")
    
    print("\n✓ All tests passed! Ready to run: shiny run src/metalotl/app.py")
    
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
