#!/usr/bin/env python3
"""
Test VTK file reading functionality for skeleton visualization.
"""
import os
import numpy as np

# Check for VTK libraries
try:
    import pyvista as pv
    PYVISTA_AVAILABLE = True
    VTK_AVAILABLE = True
    print("✅ pyvista available")
except ImportError:
    PYVISTA_AVAILABLE = False
    try:
        import vtk
        VTK_AVAILABLE = True
        print("✅ native VTK available")
    except ImportError:
        VTK_AVAILABLE = False
        print("❌ No VTK libraries available")

def test_vtk_file_reading():
    """Test reading a specific VTK file"""
    
    # Look for VTK files in the extracted_cells_conn_nonadaptive_rm directory
    base_dir = 'Y333 ATP6 ATP2/extracted_cells_conn_nonadaptive_rm'
    
    if not os.path.exists(base_dir):
        print(f"❌ Directory not found: {base_dir}")
        return
    
    # Find first available cell directory
    cell_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not cell_dirs:
        print(f"❌ No cell directories found in {base_dir}")
        return
    
    print(f"📁 Found {len(cell_dirs)} cell directories")
    
    # Test first directory
    test_dir = os.path.join(base_dir, cell_dirs[0])
    print(f"🔍 Testing directory: {test_dir}")
    
    # Look for VTK files
    vtk_files = [f for f in os.listdir(test_dir) if f.endswith('.vtk')]
    print(f"📄 Found {len(vtk_files)} VTK files: {vtk_files}")
    
    if not vtk_files:
        print("❌ No VTK files found")
        return
    
    # Test reading the first VTK file
    vtk_file = os.path.join(test_dir, vtk_files[0])
    print(f"📖 Testing file: {vtk_file}")
    
    if PYVISTA_AVAILABLE:
        try:
            mesh = pv.read(vtk_file)
            print(f"✅ pyvista read successful:")
            print(f"   Points: {mesh.n_points}")
            print(f"   Lines: {mesh.n_lines}")
            print(f"   Cells: {mesh.n_cells}")
            
            if mesh.n_lines > 0:
                print("✅ Polylines found!")
                # Try to extract polylines
                points = mesh.points
                lines = mesh.lines
                print(f"   Lines array length: {len(lines)}")
                
                # Extract first few polylines for testing
                polylines = []
                i = 0
                count = 0
                while i < len(lines) and count < 3:  # Test first 3 polylines
                    n = lines[i]
                    ids = lines[i+1:i+1+n]
                    polyline_points = points[ids]
                    polylines.append(polyline_points)
                    print(f"   Polyline {count}: {len(polyline_points)} points")
                    i += n + 1
                    count += 1
                
            else:
                print("❌ No polylines found in mesh")
                
        except Exception as e:
            print(f"❌ pyvista read failed: {e}")
    
    elif VTK_AVAILABLE:
        try:
            import vtk
            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(vtk_file)
            reader.Update()
            
            polydata = reader.GetOutput()
            points = polydata.GetPoints()
            lines = polydata.GetLines()
            
            print(f"✅ native VTK read successful:")
            print(f"   Points: {points.GetNumberOfPoints() if points else 0}")
            print(f"   Lines: {lines.GetNumberOfCells() if lines else 0}")
            
            if lines and lines.GetNumberOfCells() > 0:
                print("✅ Polylines found!")
            else:
                print("❌ No polylines found")
                
        except Exception as e:
            print(f"❌ native VTK read failed: {e}")
    
    else:
        print("❌ No VTK libraries available for testing")

if __name__ == "__main__":
    print("🧪 Testing VTK file reading functionality...")
    test_vtk_file_reading() 