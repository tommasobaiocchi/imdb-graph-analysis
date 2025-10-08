# IMDb Graph Analysis

This project implements a **social network analysis** of the IMDb dataset.  
Actors are represented as nodes in a graph, and an undirected edge is created between two actors if they have appeared in the same movie.  

The project demonstrates the use of **graph data structures and algorithms** applied to real-world data.

---

## Features

- Construction of an actor collaboration graph from IMDb datasets.
- Graph representation using adjacency lists for memory efficiency.
- Implementation of centrality measures:
  - Degree Centrality
  - Betweenness Centrality (Brandes algorithm with sampling)
  - Closeness Centrality (sampled version)
- Shortest path search between two actors (Six Degrees of Separation).
- Export of results to CSV for further analysis.

---

## Dataset

This project relies on the official [IMDb non-commercial datasets](https://developer.imdb.com/non-commercial-datasets/).  
The following files are required:

- `name.basics.tsv.gz`  
- `title.basics.tsv.gz`  
- `title.principals.tsv.gz`  

⚠️ Due to their size, these files are **not included** in this repository.  
Users should download and extract them manually before running the analysis.

---

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/tommasobaiocchi/imdb-graph-analysis.git
   ```

---

## License

This project is released under the [MIT License](LICENSE).  
It is provided for educational and research purposes only. 

