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

#### PCA plots
In the PCA plot, each individual dot represents the complete 3D landmark configuration of a specific lineage at a specific point in time.

It's like a map where similar skull shapes sit close together and different shapes sit far apart. Since we are taking timestamps every 1 million years (Ma), the plot is showing the entire history of the skull's evolution simultaneously.

- Breakdown:

-> **A Single Dot:** One lineage's skull shape at one moment in time 

-> **The Color of the Dot:** Represents when that shape existed

Yellow/Light Green dots (1 Ma) are near the root of the tree — they represent the ancestral "starting" shapes

Purple/Dark Blue dots (20 Ma) represent the "modern" shapes at the tips of the tree

The Position of the Dot: Represents what the shape looks like.

If two dots are right next to each other, those two lineages had almost identical skull shapes at those times.

If a dot is far away, that lineage evolved a very distinct, unique shape.

The Difference in Your Plots:

In the BM Plot: We can see many distinct dots spread out. This is because, in neutral evolution, lineages wander randomly. No two lineages are likely to end up in the exact same spot, so we see a wide "cloud" of individual shapes.

In the OU Plot: We see fewer distinct dots because of convergence. Selection is pulling all 30 lineages toward the exact same target shape. By the time they get to 20 Ma (purple), they have all become so similar that their dots are literally stacked on top of each other.
