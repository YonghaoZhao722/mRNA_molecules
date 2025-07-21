"""
Alpha Shape Demonstration Script

This script creates a visual demonstration of how different alpha values
affect the shape generation for irregular point patterns, helping users
understand the optimal alpha parameter for their data.
"""

import numpy as np
import matplotlib.pyplot as plt
import alphashape
from shapely.geometry import Point
from scipy.spatial import ConvexHull

def generate_sample_data(pattern='irregular'):
    """Generate sample point patterns for demonstration"""
    np.random.seed(42)
    
    if pattern == 'irregular':
        # Create an irregular biological-like pattern
        # Central cluster
        center_points = np.random.normal([0, 0], [0.3, 0.2], (30, 2))
        
        # Add some extensions/branches
        branch1 = np.random.normal([1.5, 0.5], [0.2, 0.1], (15, 2))
        branch2 = np.random.normal([-1.2, -0.8], [0.15, 0.15], (12, 2))
        branch3 = np.random.normal([0.2, 1.8], [0.1, 0.2], (10, 2))
        
        # Add some scattered points
        scattered = np.random.uniform([-2, -2], [3, 3], (8, 2))
        
        points = np.vstack([center_points, branch1, branch2, branch3, scattered])
        
    elif pattern == 'crescent':
        # Create a crescent shape
        theta = np.linspace(0, np.pi, 30)
        r1 = 2 + 0.2 * np.random.randn(30)
        r2 = 1.2 + 0.1 * np.random.randn(30)
        
        outer = np.column_stack([r1 * np.cos(theta), r1 * np.sin(theta)])
        inner = np.column_stack([r2 * np.cos(theta), r2 * np.sin(theta)])
        
        points = np.vstack([outer, inner])
        
    elif pattern == 'ring':
        # Create a ring pattern
        theta = np.linspace(0, 2*np.pi, 50)
        r = 2 + 0.3 * np.random.randn(50)
        points = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
        
        # Add some internal points
        internal = np.random.normal([0, 0], [0.5, 0.5], (10, 2))
        points = np.vstack([points, internal])
        
    return points

def plot_alpha_shape_comparison(points, alpha_values, title="Alpha Shape Comparison"):
    """Plot comparison of different alpha values"""
    n_alphas = len(alpha_values)
    fig, axes = plt.subplots(2, (n_alphas + 1) // 2, figsize=(15, 8))
    if n_alphas <= 2:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    for i, alpha in enumerate(alpha_values):
        ax = axes[i]
        
        # Plot points
        ax.scatter(points[:, 0], points[:, 1], c='blue', s=30, alpha=0.7, zorder=3)
        
        try:
            if alpha == 'convex':
                # Compute convex hull
                hull = ConvexHull(points)
                hull_points = points[hull.vertices]
                hull_points = np.vstack([hull_points, hull_points[0]])  # Close the polygon
                ax.plot(hull_points[:, 0], hull_points[:, 1], 'r-', linewidth=2, label='Convex Hull')
                ax.fill(hull_points[:, 0], hull_points[:, 1], alpha=0.2, color='red')
                ax.set_title(f'Convex Hull', fontsize=12)
            else:
                # Compute alpha shape
                alpha_shape = alphashape.alphashape(points, alpha)
                
                if alpha_shape.is_empty:
                    ax.set_title(f'α = {alpha}\n(Empty shape)', fontsize=12)
                elif hasattr(alpha_shape, 'exterior'):
                    # Single polygon
                    x, y = alpha_shape.exterior.xy
                    ax.plot(x, y, 'g-', linewidth=2, label=f'α = {alpha}')
                    ax.fill(x, y, alpha=0.2, color='green')
                    area = alpha_shape.area
                    ax.set_title(f'α = {alpha}\nArea = {area:.2f}', fontsize=12)
                elif hasattr(alpha_shape, 'geoms'):
                    # Multiple polygons
                    total_area = 0
                    for geom in alpha_shape.geoms:
                        if hasattr(geom, 'exterior'):
                            x, y = geom.exterior.xy
                            ax.plot(x, y, 'g-', linewidth=2)
                            ax.fill(x, y, alpha=0.2, color='green')
                            total_area += geom.area
                    ax.set_title(f'α = {alpha}\nTotal Area = {total_area:.2f}', fontsize=12)
                else:
                    ax.set_title(f'α = {alpha}\n(Complex shape)', fontsize=12)
                    
        except Exception as e:
            ax.set_title(f'α = {alpha}\n(Error: {str(e)[:20]})', fontsize=12)
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    # Hide unused subplots
    for i in range(n_alphas, len(axes)):
        axes[i].set_visible(False)
    
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    return fig

def demonstrate_alpha_effects():
    """Main demonstration function"""
    print("Alpha Shape Demonstration")
    print("=" * 40)
    
    # Generate different patterns
    patterns = {
        'Irregular Biological Pattern': generate_sample_data('irregular'),
        'Crescent Shape': generate_sample_data('crescent'),
        'Ring Pattern': generate_sample_data('ring')
    }
    
    # Test different alpha values
    alpha_values = [0.05, 0.1, 0.3, 0.7, 1.5, 'convex']
    
    for pattern_name, points in patterns.items():
        print(f"\nAnalyzing: {pattern_name}")
        print(f"Number of points: {len(points)}")
        
        # Create comparison plot
        fig = plot_alpha_shape_comparison(points, alpha_values, 
                                        f"{pattern_name} - Alpha Shape Analysis")
        plt.show()
        
        # Calculate optimal alpha range
        try:
            # Test a range of alpha values to find optimal
            test_alphas = np.logspace(-2, 1, 50)  # 0.01 to 10
            areas = []
            
            for alpha in test_alphas:
                try:
                    shape = alphashape.alphashape(points, alpha)
                    if hasattr(shape, 'area'):
                        areas.append(shape.area)
                    elif hasattr(shape, 'geoms'):
                        total_area = sum(geom.area for geom in shape.geoms if hasattr(geom, 'area'))
                        areas.append(total_area)
                    else:
                        areas.append(0)
                except:
                    areas.append(0)
            
            # Plot alpha vs area curve
            plt.figure(figsize=(10, 6))
            plt.subplot(1, 2, 1)
            plt.semilogx(test_alphas, areas, 'b-', linewidth=2)
            plt.xlabel('Alpha Value')
            plt.ylabel('Shape Area')
            plt.title(f'{pattern_name}\nAlpha vs Area')
            plt.grid(True, alpha=0.3)
            
            # Highlight some key alpha values
            for alpha in [0.1, 0.3, 0.7]:
                if alpha <= max(test_alphas):
                    idx = np.argmin(np.abs(test_alphas - alpha))
                    plt.plot(alpha, areas[idx], 'ro', markersize=8, 
                           label=f'α = {alpha}')
            plt.legend()
            
            # Find elbow point (optimal alpha)
            # Use second derivative to find inflection point
            if len(areas) > 10:
                areas_smooth = np.array(areas)
                # Smooth the curve
                from scipy.ndimage import gaussian_filter1d
                areas_smooth = gaussian_filter1d(areas_smooth, sigma=2)
                
                # Find second derivative
                second_deriv = np.gradient(np.gradient(areas_smooth))
                # Find point of maximum curvature
                elbow_idx = np.argmax(np.abs(second_deriv[5:-5])) + 5
                optimal_alpha = test_alphas[elbow_idx]
                
                plt.axvline(optimal_alpha, color='red', linestyle='--', 
                           label=f'Suggested α = {optimal_alpha:.3f}')
                plt.legend()
                
                print(f"Suggested optimal alpha: {optimal_alpha:.3f}")
            
            # Show detailed view of optimal alpha
            plt.subplot(1, 2, 2)
            plt.scatter(points[:, 0], points[:, 1], c='blue', s=30, alpha=0.7)
            
            try:
                if 'optimal_alpha' in locals():
                    optimal_shape = alphashape.alphashape(points, optimal_alpha)
                    if hasattr(optimal_shape, 'exterior'):
                        x, y = optimal_shape.exterior.xy
                        plt.plot(x, y, 'r-', linewidth=2, label=f'Optimal α = {optimal_alpha:.3f}')
                        plt.fill(x, y, alpha=0.2, color='red')
            except:
                pass
                
            plt.axis('equal')
            plt.title('Suggested Optimal Alpha Shape')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            print(f"Error in optimal alpha calculation: {e}")
        
        print("-" * 40)

def interactive_alpha_tuning():
    """Interactive alpha tuning tool"""
    print("\nInteractive Alpha Tuning")
    print("=" * 30)
    
    # Generate sample data
    points = generate_sample_data('irregular')
    
    while True:
        try:
            alpha_input = input("\nEnter alpha value (or 'quit' to exit): ")
            if alpha_input.lower() == 'quit':
                break
                
            alpha = float(alpha_input)
            
            # Compute and display alpha shape
            plt.figure(figsize=(8, 6))
            plt.scatter(points[:, 0], points[:, 1], c='blue', s=30, alpha=0.7, label='Data Points')
            
            try:
                alpha_shape = alphashape.alphashape(points, alpha)
                
                if alpha_shape.is_empty:
                    print(f"Alpha {alpha} resulted in an empty shape. Try a larger value.")
                elif hasattr(alpha_shape, 'exterior'):
                    x, y = alpha_shape.exterior.xy
                    plt.plot(x, y, 'g-', linewidth=2, label=f'Alpha Shape (α = {alpha})')
                    plt.fill(x, y, alpha=0.2, color='green')
                    print(f"Alpha {alpha}: Area = {alpha_shape.area:.3f}")
                elif hasattr(alpha_shape, 'geoms'):
                    total_area = 0
                    for i, geom in enumerate(alpha_shape.geoms):
                        if hasattr(geom, 'exterior'):
                            x, y = geom.exterior.xy
                            plt.plot(x, y, 'g-', linewidth=2, 
                                   label=f'Component {i+1}' if i < 3 else '')
                            plt.fill(x, y, alpha=0.2, color='green')
                            total_area += geom.area
                    print(f"Alpha {alpha}: {len(alpha_shape.geoms)} components, Total area = {total_area:.3f}")
                
            except Exception as e:
                print(f"Error computing alpha shape: {e}")
            
            plt.axis('equal')
            plt.title(f'Alpha Shape with α = {alpha}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()
            
        except ValueError:
            print("Please enter a valid number or 'quit'")
        except KeyboardInterrupt:
            break

if __name__ == "__main__":
    # Check if required packages are available
    try:
        import alphashape
        import scipy.ndimage
    except ImportError as e:
        print(f"Missing required package: {e}")
        print("Please install: pip install alphashape scipy")
        exit(1)
    
    # Run demonstrations
    print("Welcome to Alpha Shape Demonstration!")
    print("\nThis script will show you how different alpha values affect shape generation.")
    print("Close each plot window to continue to the next demonstration.")
    
    # Main demonstration
    demonstrate_alpha_effects()
    
    # Interactive tuning
    response = input("\nWould you like to try interactive alpha tuning? (y/n): ")
    if response.lower().startswith('y'):
        interactive_alpha_tuning()
    
    print("\nDemo complete! Use these insights to choose optimal alpha values for your data.") 