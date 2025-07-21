#!/usr/bin/env python3
"""
Check if pyvista is available for VTK skeleton visualization.
"""

try:
    import pyvista as pv
    print("✅ pyvista is available!")
    print(f"   Version: {pv.__version__}")
    
    # Test basic functionality
    try:
        # Create a simple test mesh
        mesh = pv.Sphere()
        print(f"   Basic functionality test: OK ({mesh.n_points} points)")
    except Exception as e:
        print(f"   Warning: Basic functionality test failed: {e}")
        
except ImportError:
    print("❌ pyvista is NOT available!")
    print("\nTo install pyvista, run:")
    print("   pip install pyvista")
    print("\nOr if using conda:")
    print("   conda install -c conda-forge pyvista")
    print("\nNote: pyvista is required for proper skeleton visualization")
    print("      (reading VTK files with correct topology)")

print("\nTesting complete.") 