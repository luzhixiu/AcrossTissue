import matplotlib.pyplot as plt

import seaborn as sns; sns.set(color_codes=True)
iris = sns.load_dataset("iris")
species = iris.pop("species")

print(type(iris))
g = sns.clustermap(iris)
plt.show()