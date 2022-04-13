import open3d as o3d
import numpy as np

# xyz = np.random.rand(100, 3)
# pcd = o3d.geometry.PointCloud()
# pcd.points = o3d.utility.Vector3dVector(xyz)
# o3d.io.write_point_cloud("./data.ply", pcd)
# o3d.visualization.draw_geometries([pcd])

pcd = o3d.io.read_point_cloud('./newPC.txt', format='xyzrgb')
sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0, resolution=20)
o3d.io.write_point_cloud("./data.ply", pcd)
print(len(pcd.points))
o3d.visualization.draw_geometries([pcd, sphere])

# o3d.visualization.draw({'name': 'test', 'geometry': [sphere], 'material': mat})