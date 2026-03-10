# Morphometric-traits-over-time

1. Why use landmarks ?

Landmarks allow to capture the spatial relationship between parts of the bone. When no DNA data are available, which is common in paleontology or archaeology, landmarks act as a "proxy" for biological relatedness. While not perfect, the shape of a bone is a physical expression of an organism's genes. By comparing landmark configurations, we can build "shape trees" (phylogenies) that often mirror genetic trees.

Once all landmarks are required, they are turned into a Principal Component Analysis (PCA). This creates a visual "map" where each dot is a different specimen. If the dots cluster together, their bone shapes are nearly identical, suggesting they are closely related or share the same environment.

Pipeline: simulate skull 3D shape by simulating landmark configurations. Each shape is a set of *k* landmarks in *p* directions (here 3 directions). 

Stages:
1. Create the 3D cranial landmark configuration of the ancestor
2. Use the `pbtree` function of the `phytools` package, to generate a tree that will represent the time scale
3. Simulate the evolution of the landmark coordinates across the branches of the phylogenetic tree using the `mvSIM` function of the `mvMORPH` package. This simulation of the traits should be applied with 2 different models: 
- **Neutral evolution** = Brownian Motion, where model traits evolve randomly over time
- **Selection** = Ornstein–Uhlenbeck, where the traits evolve towards an adaptive optimum.
4. Output skull configurations at every 1 million years using the `fastAnc` function in `phytools`
5. Compare the "shape distances" to the "evolutionary time distances"
