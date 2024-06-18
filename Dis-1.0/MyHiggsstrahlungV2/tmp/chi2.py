import yoda
import matplotlib.pyplot as plt

h1 = yoda.read("Analysisbgwzc.yoda", asdict=False)

plt.hist(h1)

plt.show