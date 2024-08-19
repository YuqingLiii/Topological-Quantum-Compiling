import numpy as np
from timeit import default_timer as timer
def Rz(angle):
    """Rotation around the Z-axis."""
    return np.array([
        [np.exp(-0.5j * angle), 0],
        [0, np.exp(0.5j * angle)]
    ], dtype=complex)

def gate_distance(guess_coeff, target_matrix):
    a, b, c, d = guess_coeff
    tau = (np.sqrt(5) - 1) / 2
    F = np.array([
        [tau, np.sqrt(tau)],
        [np.sqrt(tau), -tau]
    ], dtype=complex)
    
    estim_matrix = np.exp(1j * a) * Rz(b) @ F @ Rz(c) @ F @ Rz(d)
    
    difference_matrix = estim_matrix - target_matrix
    frobenius_norm = np.linalg.norm(difference_matrix, 'fro')
    
    return frobenius_norm

def gradient_descent(target_matrix, learning_rate=0.0001, max_iter=100000):
    # Initialize coefficients (start with random values)
    guess_coeff = np.random.rand(4) * 2 * np.pi
    
    for i in range(max_iter):
        # Compute the gradient numerically
        gradients = np.zeros_like(guess_coeff, dtype=float)
        epsilon = 1e-8
        
        for j in range(len(guess_coeff)):
            coeffs_step_forward = np.copy(guess_coeff)
            coeffs_step_forward[j] += epsilon
            loss_forward = gate_distance(coeffs_step_forward, target_matrix)
            
            coeffs_step_backward = np.copy(guess_coeff)
            coeffs_step_backward[j] -= epsilon
            loss_backward = gate_distance(coeffs_step_backward, target_matrix)
            
            gradients[j] = (loss_forward - loss_backward) / (2 * epsilon)
        
        # Update the coefficients using gradient descent
        guess_coeff -= learning_rate * gradients
        
        # Calculate the current loss (distance)
        current_loss = gate_distance(guess_coeff, target_matrix)
        
        # Print the progress
        if i % 100 == 0:
            print(f"Iteration {i}: Loss = {current_loss}")
        
        # Convergence check (optional)
        if current_loss < 1e-6:
            print("Converged.")
            break
    
    return guess_coeff

if __name__ == '__main__':
    time1=timer()
    # Define the target matrix (e.g., Hadamard matrix)
    target_matrix = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # Run the gradient descent algorithm
    optimal_coeffs = gradient_descent(target_matrix)
    time2=timer()
    print(f"Optimal coefficients: {optimal_coeffs}, using {time2-time1} s")
