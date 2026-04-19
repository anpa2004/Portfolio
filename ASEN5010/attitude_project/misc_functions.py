import numpy as np

def vector_extract(vec_list:list[float])->list[float]:
    """ 
    Take a list of 3x1 vectors and return lists of components

    Args:
        vec_list: list of vectors to extract components

    Returns:
        a,b,c: components to return as lists
    """

    a = [v[0] for v in vec_list]
    b = [v[1] for v in vec_list]
    c = [v[2] for v in vec_list]
    return a,b,c

def create_submission_txt(filename: str, value) -> None:
    """ 
    Create a submission text file with forced decimal formatting (no scientific notation).
    """
    
    # Helper to format any number into a long decimal string
    # Change '.12f' to '.15f' if you need even more precision
    def fmt(x): 
        return f"{x:.12f}"

    with open(filename, "w") as f:
        # Case 1: Scalar (float or int)
        if isinstance(value, (float, int, np.number)):
            f.write(fmt(value) + "\n")
            
        # Case 2: List
        elif isinstance(value, list):
            formatted_list = [fmt(x) for x in value]
            f.write(" ".join(formatted_list) + "\n")
            
        # Case 3: NumPy Array
        elif isinstance(value, np.ndarray):
            # Flatten and format every element individually
            formatted_vals = [fmt(x) for x in value.flatten()]
            f.write(" ".join(formatted_vals) + "\n")

def unitize(vec:list[float])->list[float]:
    """ 
    Take a unit vector and return the unit vector in the same direction
    """
    return vec/np.linalg.norm(vec)

def print_latex(a:list[float]) -> None:
    """ 
    Print matrix as a latex matrix for simplicity (3x3) or (3,1)
    """
    if a.shape == (3,3):
        print('\\begin{bmatrix}')
        print(f'{a[0,0]} & {a[0,1]} & {a[0,2]}\\\\')
        print(f'{a[1,0]} & {a[1,1]} & {a[1,2]}\\\\')
        print(f'{a[2,0]} & {a[2,1]} & {a[2,2]}\\\\')
        print('\\end{bmatrix}')
    elif a.shape == (3,1):
        print('\\begin{bmatrix}')
        print(f'{a[0][0]} \\\\ {a[1][0]} \\\\ {a[2][0]}')
        print('\\end{bmatrix}')
    elif a.shape == (1,3) or a.shape == (3,):
        print('\\begin{bmatrix}')
        print(f'{a[0]} \\\\ {a[1]} \\\\ {a[2]}')
        print('\\end{bmatrix}')
    else: 
        raise TypeError(f'Incompatible shape for this function {a.shape}')
    
def un_tilde(v_tilde:list)->list:
        """   
        Take a tilde matrix and return the corresponding vector
        """
        v = np.array([v_tilde[2,1],v_tilde[0,2],v_tilde[1,0]])
        return v

def tilde(w):
    return np.array([
        [0, -w[2], w[1]],
        [w[2], 0, -w[0]],
        [-w[1], w[0], 0]
    ])