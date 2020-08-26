from Bio.Phylo.TreeConstruction import DistanceCalculator,ParsimonyScorer,NNITreeSearcher,ParsimonyTreeConstructor, DistanceTreeConstructor
from Bio import AlignIO, Phylo

# Read from file and calculate distance matrics
alignment = AlignIO.read('example.txt', 'fasta')
calculator = DistanceCalculator('identity')
distance_matrics = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()

nj_tree = constructor.nj(distance_matrics)
upgma_tree = constructor.upgma(distance_matrics)

print(nj_tree)
Phylo.draw(nj_tree)

print(upgma_tree)
Phylo.draw(upgma_tree)


score_pars_tree = ParsimonyScorer()
searcher_pars_tree = NNITreeSearcher(score_pars_tree)
constructor_pars_tree = ParsimonyTreeConstructor(searcher_pars_tree, nj_tree)
pars_tree = constructor_pars_tree.build_tree(alignment)

print(pars_tree)
Phylo.draw(pars_tree)



