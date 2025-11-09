import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R
import os

def main(scene_file : str, use_gui : bool = True):
	# Set up the simulator
	base = sph.Exec.SimulatorBase()
	base.init(useGui=use_gui, sceneFile=scene_file)
	
	sim = sph.Simulation.getCurrent()

	# Create an imgui simulator
	gui = sph.GUI.Simulator_GUI_imgui(base)
	base.setGui(gui)
	
	base.initSimulation()
		
	base.runSimulation()
	base.cleanup()

if __name__ == "__main__":
    import meshio
    import argparse

    parser = argparse.ArgumentParser("pysplishsplash")
    parser.add_argument("--scene-file", type=str, default="data/Scenes/fix_right_0p5m_ym5e6_sph.json")
    parser.add_argument("--mesh-file", type=str, default="data/models/beam5m_res50x10x10.msh")
    parser.add_argument("--use-gui", type=bool, default=True)
    args = parser.parse_args()

    scene_file = args.scene_file
    mesh_file = args.mesh_file
    use_gui = args.use_gui
    # Load the mesh data from the PLY file
    mesh = meshio.read(mesh_file)

    # Access the vertex coordinates
    vertices = mesh.points
    print(f"Vertices shape: {vertices.shape}")
    print(f"Mesh: {mesh}")

    particle_filename = mesh_file.replace(".msh", "_particles.bgeo")
    if (os.path.exists(particle_filename)):
        print(f"Particle file already exists: {particle_filename}. Overwriting...")
    positions_np = np.asarray(mesh.points, dtype=np.float32)
    # positions = [positions_np.tolist() for row in positions_np]
    sph.Utilities.PartioReaderWriter.writeParticles(
        particle_filename,
        positions_np,
        None,
        0.0,
    )

    main(scene_file, use_gui)
