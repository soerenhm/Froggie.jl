# Uniform vector
A = [1, 2, 3, 4]
@test binalong(A, 2) == [6., 4.]
x = [0, 2, 4, 8]
@test binalong(A, x, 3) == [3., 3., 4.]
B = [1, 1]
@test binalong(B, 4) ≈ [1.0/3, 1.0/3, 1.0/3, 1.]
# Matrix
C = [1 1 1; 2 2 2; 3 3 3]
@test binalong(C, 2, 2) == [2. 1; 4. 2; 6. 3]
@test binalong(C, 1, 2) == [3. 3. 3; 3. 3 3]
@test binalong(C', 1, 2) == binalong(C, 2, 2)'
# AxisArray
D = AxisArray(C; x=0:2, y=1:3)