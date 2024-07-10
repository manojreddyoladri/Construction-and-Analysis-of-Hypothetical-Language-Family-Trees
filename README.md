# Construction-and-Analysis-of-Hypothetical-Language-Family-Trees

# Description

This project implements two algorithms, UPGMA (Unweighted Pair Group Method with Arithmetic Mean) and CMM (Character-based Matrix Method), to build hypothetical language family trees using the Hw1_words data file. The project aims to compare the results of these two methods and analyze their effectiveness in explaining and reflecting the actual history of language development.

# Installation

1. Clone the repository:
```
git clone https://github.com/manojreddyoladri/Construction-and-Analysis-of-Hypothetical-Language-Family-Trees.git
cd Construction-and-Analysis-of-Hypothetical-Language-Family-Trees
```
2. Install the required dependencies:
```
pip install -r requirements.txt
```
# Usage
1. To run the UPGMA and CMM algorithm:
```
python main.py
```
2. The resulting language family trees will be saved as output files in the Plots directory.

# Features
1. UPGMA Algorithm: Implements the UPGMA algorithm to build a hypothetical language family tree.
2. CMM Algorithm: Implements the CMM algorithm to build a hypothetical language family tree.

# Technologies Used
1. Python
2. Openpyxl
3. Numpy
4. Biopython
5. Matplotlib

# Analysis
1. Detailed Explanation: The CMM algorithm provides a more detailed explanation of the development of the language family tree. The CMM algorithm considers the shared mutations between languages to construct the family tree. Using a hierarchical clustering method, the algorithm provides a very detailed phylogenetic tree, which displays the sequence of mutations as well as the relationships between languages. Essentially, it can show us the order for which these mutations occurred; something that the UPGMA algorithm cannot. For matrices with rare mutations, it can be more efficient to count those mutation occurrences rather than calculating their distances, and the CMM algorithm is especially encouraged in these situations.
2. Accuracy: The CMM provides more accurate family tree, particularly in datasets characterized by long branch lengths, such as those seen in Hungarian and Mansi. Long branch lengths are cases where languages diverged a long time ago – though they may still have substantial commonalities. The UPGMA algorithm will only consider the extent to which two languages are similar, and may not be well suited to distinguishing languages that diverted recently from those that diverted long ago. In contrast, CMM considers the extent to which two languages have common differences (or mutations). Common mutations suggest a common history of divergence between languages. Intuitively, languages that diverted from a common ancestor at the same time (only diverting themselves later) should share common mutations indicative of that initial diversion. For instance, the UPGMA closely relates Saami to Hungarian and Mansi. This relationship does not appear in the CMM analysis because, although Saami has many commonalities with Hungarian and Mansi, it has few common mutations. CMM takes into account common mutations, a more precise depiction of linguistic relationships in scenarios with substantial variation in the data.
