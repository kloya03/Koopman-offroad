function [normalized_closeness,closeness] = MatrixSimilarity(A,B)

difference = A - B; 
closeness = norm(difference, 'fro'); 
normalized_closeness = closeness / norm(Xi_N1, 'fro'); 