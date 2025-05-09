# Network Lasso with Splitting Methods: STAT 8054 Final Project

This project implements and compares two splitting methods—ADMM and AMA—for solving the Network Lasso problem on synthetic data generated from a Tau protein model. The goal is to evaluate their convergence behavior, computational efficiency, and clustering performance under different regularization parameters.

## Project Structure

| File/Folder                 | Description                                                     |
|-----------------------------|-----------------------------------------------------------------|
| `8054_project_slides.pdf`   | Presentation slides                                             |
| `8054_project_write_up.pdf` | Final report                                                    |
| `experiments/`              | R scripts for simulation experiments                            |
| `results/`                  | R scripts to load summary dataframes and generate visualizations|
| `README.md`                 | This project description                                        |
| `Network_Lasso.pdf`         | reference 1                                                     |
| `Splitting_Methods.pdf`     | reference 2                                                     |

## How to Reproduce

1. Run `task.R` in each subfolder of `experiments/`. This will produce two CSV files for ADMM and AMA, respectively.
2. Run the R scripts in the `results/` folder to generate the corresponding plots.

## Main Results

- AMA converges faster during the early runtime for both small and large \(\lambda\) values but struggles during the merging phase.
- ADMM maintains stable convergence across the full regularization path.
- Warm start accelerates both methods but does not fully resolve AMA’s convergence issues.

## References

- Hallac et al. (2015). *Network Lasso: Clustering and Optimization in Large Graphs.*
- Chi & Lange (2015). *Splitting Methods for Convex Clustering.*
