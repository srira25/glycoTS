# mi_computation.py
from sklearn.feature_selection import mutual_info_regression
from scipy.stats import pearsonr, spearmanr

def compute_mi(gene_expression, tf_gene_expression):
    
    tf_gene_expression = tf_gene_expression.reshape(-1, 1)
    mi = mutual_info_regression(tf_gene_expression, gene_expression)[0]
    
    return mi  

def compute_pearson(gene_expression, tf_gene_expression):
    
    pearson_corr, _ = pearsonr(tf_gene_expression.ravel(), gene_expression.ravel())
    return pearson_corr

def compute_spearman(gene_expression, tf_gene_expression):
    
    spearman_corr, _ = spearmanr(tf_gene_expression.ravel(), gene_expression.ravel())
    return spearman_corr