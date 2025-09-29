#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <algorithm>
#include <random>
#include <chrono>
#include <functional>

using namespace std;

class IMDBGraph {
private:
    // Internal ID mapping for memory optimization
    unordered_map<string, int> actor_id_to_int;
    unordered_map<string, int> movie_id_to_int;
    vector<string> int_to_actor_id;
    vector<string> int_to_movie_id;
    
    struct Actor {
        int id;
        string name;
        vector<int> movie_ids;
        vector<int> collaborators; // Changed from unordered_set to vector for memory efficiency
        
        Actor() : id(-1) {}
        Actor(int actor_id, const string& actor_name) 
            : id(actor_id), name(actor_name) {}
    };
    
    struct Movie {
        int id;
        string title;
        int year;
        vector<int> actor_ids;
        
        Movie() : id(-1), year(0) {}
        Movie(int movie_id, const string& movie_title, int movie_year)
            : id(movie_id), title(movie_title), year(movie_year) {}
    };

    vector<Actor> actors;
    vector<Movie> movies;
    
    // Internal ID management
    int getOrCreateActorID(const string& actor_str_id) {
        auto it = actor_id_to_int.find(actor_str_id);
        if (it != actor_id_to_int.end()) {
            return it->second;
        }
        int new_id = int_to_actor_id.size();
        actor_id_to_int[actor_str_id] = new_id;
        int_to_actor_id.push_back(actor_str_id);
        actors.emplace_back(new_id, "");
        return new_id;
    }
    
    int getOrCreateMovieID(const string& movie_str_id) {
        auto it = movie_id_to_int.find(movie_str_id);
        if (it != movie_id_to_int.end()) {
            return it->second;
        }
        int new_id = int_to_movie_id.size();
        movie_id_to_int[movie_str_id] = new_id;
        int_to_movie_id.push_back(movie_str_id);
        movies.emplace_back(new_id, "", 0);
        return new_id;
    }
    
    // Check if actor has acting profession
    bool isActor(const string& professions) {
        if (professions == "\\N") return false;
        
        vector<string> tokens;
        stringstream ss(professions);
        string token;
        while (getline(ss, token, ',')) {
            if (token == "actor" || token == "actress") {
                return true;
            }
        }
        return false;
    }

public:
    IMDBGraph() {
        // Enable fast I/O
        ios::sync_with_stdio(false);
        cin.tie(nullptr);
    }
    
    // Parse TSV line into fields
    vector<string> splitTSVLine(const string& line) {
        vector<string> fields;
        stringstream ss(line);
        string field;
        
        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }
        return fields;
    }
    
    // Parse year string, handle missing values
    int parseYear(const string& year_str) {
        if (year_str == "\\N" || year_str.empty()) {
            return -1;
        }
        try {
            return stoi(year_str);
        } catch (const exception& e) {
            return -1;
        }
    }
    
    // Load actors from name.basics.tsv
    bool loadActorsFromTSV(const string& filename) {
        cout << "Loading actors from: " << filename << endl;
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open file " << filename << endl;
            return false;
        }
        
        string line;
        getline(file, line); // Skip header
        
        int count = 0;
        while (getline(file, line)) {
            auto fields = splitTSVLine(line);
            if (fields.size() < 6) continue;
            
            string actor_str_id = fields[0];
            string name = fields[1];
            string professions = fields[4];
            
            // Filter only actors/actresses using exact profession matching
            if (isActor(professions)) {
                int actor_id = getOrCreateActorID(actor_str_id);
                actors[actor_id].name = name;
                
                // Parse knownForTitles
                string movie_list = fields[5];
                if (movie_list != "\\N") {
                    stringstream movie_ss(movie_list);
                    string movie_str_id;
                    while (getline(movie_ss, movie_str_id, ',')) {
                        actors[actor_id].movie_ids.push_back(getOrCreateMovieID(movie_str_id));
                    }
                }
                count++;
            }
            
            if (count % 10000 == 0) {
                cout << "Loaded " << count << " actors..." << endl;
            }
        }
        
        cout << "Total actors loaded: " << count << endl;
        return true;
    }
    
    // Load movies from title.basics.tsv
    bool loadMoviesFromTSV(const string& filename) {
        cout << "Loading movies from: " << filename << endl;
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open file " << filename << endl;
            return false;
        }
        
        string line;
        getline(file, line); // Skip header
        
        int count = 0;
        while (getline(file, line)) {
            auto fields = splitTSVLine(line);
            if (fields.size() < 9) continue;
            
            string movie_str_id = fields[0];
            string type = fields[1];
            string title = fields[2];
            int year = parseYear(fields[5]);
            
            // Consider only movies (not TV series) with valid year
            if (type == "movie" && year > 1900) {
                int movie_id = getOrCreateMovieID(movie_str_id);
                movies[movie_id].title = title;
                movies[movie_id].year = year;
                count++;
            }
            
            if (count % 10000 == 0) {
                cout << "Loaded " << count << " movies..." << endl;
            }
        }
        
        cout << "Total movies loaded: " << count << endl;
        return true;
    }
    
    // Load actor-movie relationships from title.principals.tsv
    bool loadPrincipalsFromTSV(const string& filename) {
        cout << "Loading principals from: " << filename << endl;
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open file " << filename << endl;
            return false;
        }
        
        string line;
        getline(file, line); // Skip header
        
        int count = 0;
        while (getline(file, line)) {
            auto fields = splitTSVLine(line);
            if (fields.size() < 6) continue;
            
            string movie_str_id = fields[0];
            string actor_str_id = fields[2];
            string category = fields[3];
            
            // Consider only actors/actresses in valid movies
            if ((category == "actor" || category == "actress")) {
                auto movie_it = movie_id_to_int.find(movie_str_id);
                auto actor_it = actor_id_to_int.find(actor_str_id);
                
                if (movie_it != movie_id_to_int.end() && actor_it != actor_id_to_int.end()) {
                    int movie_id = movie_it->second;
                    int actor_id = actor_it->second;
                    
                    if (movie_id < movies.size() && actor_id < actors.size()) {
                        movies[movie_id].actor_ids.push_back(actor_id);
                        actors[actor_id].movie_ids.push_back(movie_id);
                        count++;
                    }
                }
            }
            
            if (count % 100000 == 0) {
                cout << "Processed " << count << " roles..." << endl;
            }
        }
        
        cout << "Total roles loaded: " << count << endl;
        return true;
    }
    
    // Build collaboration graph where actors are connected if they worked together
    void buildCollaborationGraph() {
        cout << "Building collaboration graph..." << endl;
        auto start = chrono::high_resolution_clock::now();
        
        int collaborations_created = 0;
        
        // For each movie, connect all actors who worked together
        for (const auto& movie : movies) {
            if (movie.id == -1) continue;
            
            const auto& actor_ids = movie.actor_ids;
            
            // Create complete clique between all actors in the movie
            for (size_t i = 0; i < actor_ids.size(); ++i) {
                for (size_t j = i + 1; j < actor_ids.size(); ++j) {
                    int actor1 = actor_ids[i];
                    int actor2 = actor_ids[j];
                    
                    if (actor1 < actors.size() && actor2 < actors.size()) {
                        // Add bidirectional collaboration
                        actors[actor1].collaborators.push_back(actor2);
                        actors[actor2].collaborators.push_back(actor1);
                        collaborations_created++;
                    }
                }
            }
        }
        
        // Remove duplicates and optimize memory
        cout << "Optimizing memory by removing duplicates..." << endl;
        for (auto& actor : actors) {
            if (actor.id == -1) continue;
            
            // Sort and remove duplicates
            sort(actor.collaborators.begin(), actor.collaborators.end());
            auto last = unique(actor.collaborators.begin(), actor.collaborators.end());
            actor.collaborators.erase(last, actor.collaborators.end());
            
            // Shrink to fit
            vector<int>(actor.collaborators).swap(actor.collaborators);
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        
        cout << "Graph built in " << duration.count() << "ms" << endl;
        cout << "Total collaborations: " << collaborations_created << endl;
        cout << "Unique edges: " << getNumEdges() << endl;
    }
    
    // BFS to find shortest path between two actors (Six Degrees of Separation)
    vector<int> findShortestPath(int actor1_id, int actor2_id) {
        if (actor1_id == actor2_id) {
            return {actor1_id};
        }
        
        auto start = chrono::high_resolution_clock::now();
        
        // BFS initialization
        vector<int> parent(actors.size(), -1);
        queue<int> q;
        vector<bool> visited(actors.size(), false);
        
        q.push(actor1_id);
        visited[actor1_id] = true;
        parent[actor1_id] = -1;
        
        bool found = false;
        
        while (!q.empty() && !found) {
            int current = q.front();
            q.pop();
            
            // Explore all collaborators
            for (int collaborator : actors[current].collaborators) {
                if (!visited[collaborator]) {
                    visited[collaborator] = true;
                    parent[collaborator] = current;
                    q.push(collaborator);
                    
                    // Check if we reached destination
                    if (collaborator == actor2_id) {
                        found = true;
                        break;
                    }
                }
            }
            if (found) break;
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        
        if (!found) {
            cout << "No path found between the actors" << endl;
            return {};
        }
        
        // Reconstruct path
        vector<int> path;
        int node = actor2_id;
        while (node != -1) {
            path.push_back(node);
            node = parent[node];
        }
        reverse(path.begin(), path.end());
        
        cout << "Path found in " << duration.count() << "ms" << endl;
        cout << "Path length: " << path.size() - 1 << " degrees" << endl;
        
        return path;
    }
    
    // Degree Centrality: Number of direct connections
    vector<int> computeDegreeCentrality() {
        cout << "Computing Degree Centrality..." << endl;
        auto start = chrono::high_resolution_clock::now();
        
        vector<int> centrality(actors.size(), 0);
        
        for (int i = 0; i < actors.size(); ++i) {
            if (actors[i].id != -1) {
                centrality[i] = actors[i].collaborators.size();
            }
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Degree Centrality computed in " << duration.count() << "ms" << endl;
        
        return centrality;
    }
    
    // Betweenness Centrality using Brandes algorithm with sampling
    vector<double> computeBetweennessCentrality(int sample_size = 1000) {
        cout << "Computing Betweenness Centrality (Brandes algorithm with sampling: " << sample_size << ")..." << endl;
        cout << "Note: Values are averaged over samples for relative ranking" << endl;
        auto start = chrono::high_resolution_clock::now();
        
        vector<double> centrality(actors.size(), 0.0);
        
        // Get valid actor indices
        vector<int> valid_actors;
        for (int i = 0; i < actors.size(); ++i) {
            if (actors[i].id != -1 && !actors[i].collaborators.empty()) {
                valid_actors.push_back(i);
            }
        }
        
        if (valid_actors.empty()) return centrality;
        
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, valid_actors.size() - 1);
        
        for (int sample = 0; sample < sample_size; ++sample) {
            int s = valid_actors[dis(gen)];
            
            stack<int> S;
            vector<vector<int>> P(actors.size());
            vector<int> sigma(actors.size(), 0);
            vector<int> dist(actors.size(), -1);
            
            sigma[s] = 1;
            dist[s] = 0;
            queue<int> Q;
            Q.push(s);
            
            while (!Q.empty()) {
                int v = Q.front(); Q.pop();
                S.push(v);
                
                for (int w : actors[v].collaborators) {
                    if (dist[w] < 0) {
                        dist[w] = dist[v] + 1;
                        Q.push(w);
                    }
                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        P[w].push_back(v);
                    }
                }
            }
            
            vector<double> delta(actors.size(), 0.0);
            while (!S.empty()) {
                int w = S.top(); S.pop();
                for (int v : P[w]) {
                    delta[v] += (sigma[v] * 1.0 / sigma[w]) * (1.0 + delta[w]);
                }
                if (w != s) {
                    centrality[w] += delta[w];
                }
            }
            
            if ((sample + 1) % 100 == 0) {
                cout << "Processed " << (sample + 1) << " samples..." << endl;
            }
        }
        
        // Average over samples (not normalized for undirected graph)
        for (int i = 0; i < centrality.size(); ++i) {
            centrality[i] /= (sample_size * 1.0);
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Betweenness Centrality computed in " << duration.count() << "ms" << endl;
        
        return centrality;
    }
    
    // Closeness Centrality: Average distance to all other actors (sampled version)
    vector<double> computeClosenessCentrality(int sample_size = 1000) {
        cout << "Computing Closeness Centrality (sampled version, " << sample_size << " sources)..." << endl;
        cout << "Note: Only sampled sources have closeness > 0" << endl;
        auto start = chrono::high_resolution_clock::now();
        
        vector<double> closeness(actors.size(), 0.0);
        
        // Get valid actor indices
        vector<int> valid_actors;
        for (int i = 0; i < actors.size(); ++i) {
            if (actors[i].id != -1 && !actors[i].collaborators.empty()) {
                valid_actors.push_back(i);
            }
        }
        
        if (valid_actors.empty()) return closeness;
        
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, valid_actors.size() - 1);
        
        for (int i = 0; i < sample_size; ++i) {
            int source = valid_actors[dis(gen)];
            
            vector<int> dist(actors.size(), -1);
            queue<int> q;
            dist[source] = 0;
            q.push(source);
            
            double total_distance = 0.0;
            int reachable_count = 0;
            
            while (!q.empty()) {
                int current = q.front(); q.pop();
                
                for (int neighbor : actors[current].collaborators) {
                    if (dist[neighbor] == -1) {
                        dist[neighbor] = dist[current] + 1;
                        total_distance += dist[neighbor];
                        reachable_count++;
                        q.push(neighbor);
                    }
                }
            }
            
            if (reachable_count > 0) {
                closeness[source] = reachable_count / total_distance;
            }
            
            if ((i + 1) % 100 == 0) {
                cout << "Processed " << (i + 1) << " samples..." << endl;
            }
        }
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Closeness Centrality computed in " << duration.count() << "ms" << endl;
        
        return closeness;
    }
    
    // Get top N actors by centrality measure
    vector<pair<int, double>> getTopCentrality(const vector<double>& centrality, int top_n = 10) {
        vector<pair<int, double>> sorted_centrality;
        
        for (int i = 0; i < centrality.size(); ++i) {
            if (actors[i].id != -1 && centrality[i] > 0) {
                sorted_centrality.emplace_back(i, centrality[i]);
            }
        }
        
        sort(sorted_centrality.begin(), sorted_centrality.end(),
            [](const pair<int, double>& a, const pair<int, double>& b) {
                return a.second > b.second;
            });
        
        if (top_n > 0 && top_n < sorted_centrality.size()) {
            sorted_centrality.resize(top_n);
        }
        
        return sorted_centrality;
    }
    
    vector<pair<int, int>> getTopDegree(const vector<int>& degree_centrality, int top_n = 10) {
        vector<pair<int, int>> sorted_degree;
        
        for (int i = 0; i < degree_centrality.size(); ++i) {
            if (actors[i].id != -1 && degree_centrality[i] > 0) {
                sorted_degree.emplace_back(i, degree_centrality[i]);
            }
        }
        
        sort(sorted_degree.begin(), sorted_degree.end(),
            [](const pair<int, int>& a, const pair<int, int>& b) {
                return a.second > b.second;
            });
        
        if (top_n > 0 && top_n < sorted_degree.size()) {
            sorted_degree.resize(top_n);
        }
        
        return sorted_degree;
    }
    
    // Save results to CSV file
    void saveResultsToCSV(const vector<pair<int, double>>& results, const string& filename, const string& centrality_type) {
        ofstream file(filename);
        file << "ActorName," << centrality_type << "Score\n";
        
        for (const auto& [actor_id, score] : results) {
            file << "\"" << actors[actor_id].name << "\"," << score << "\n";
        }
        
        file.close();
        cout << "Results saved to " << filename << endl;
    }
    
    void saveDegreeResultsToCSV(const vector<pair<int, int>>& results, const string& filename) {
        ofstream file(filename);
        file << "ActorName,DegreeCentrality\n";
        
        for (const auto& [actor_id, degree] : results) {
            file << "\"" << actors[actor_id].name << "\"," << degree << "\n";
        }
        
        file.close();
        cout << "Degree results saved to " << filename << endl;
    }
    
    // Utility functions
    size_t getNumEdges() const {
        size_t total_edges = 0;
        for (const auto& actor : actors) {
            if (actor.id != -1) {
                total_edges += actor.collaborators.size();
            }
        }
        return total_edges / 2; // Each edge counted twice
    }
    
    string getActorName(int actor_id) const {
        if (actor_id < 0 || actor_id >= actors.size() || actors[actor_id].id == -1) {
            return "Unknown";
        }
        return actors[actor_id].name;
    }
    
    int findActorByName(const string& name) {
        for (int i = 0; i < actors.size(); ++i) {
            if (actors[i].id != -1 && actors[i].name.find(name) != string::npos) {
                return i;
            }
        }
        return -1;
    }
    
    int getActorID(const string& actor_str_id) const {
        auto it = actor_id_to_int.find(actor_str_id);
        return it != actor_id_to_int.end() ? it->second : -1;
    }
    
    void printGraphStats() const {
        cout << "\n=== IMDB GRAPH STATISTICS ===" << endl;
        cout << "Actors: " << countValidActors() << endl;
        cout << "Movies: " << countValidMovies() << endl;
        cout << "Edges (collaborations): " << getNumEdges() << endl;
        
        // Calculate average degree
        double avg_degree = 0;
        int max_degree = 0;
        string most_connected_actor;
        int valid_actors = 0;
        
        for (const auto& actor : actors) {
            if (actor.id != -1) {
                int degree = actor.collaborators.size();
                avg_degree += degree;
                valid_actors++;
                if (degree > max_degree) {
                    max_degree = degree;
                    most_connected_actor = actor.name;
                }
            }
        }
        
        if (valid_actors > 0) {
            avg_degree /= valid_actors;
            cout << "Average degree: " << avg_degree << endl;
            cout << "Maximum degree: " << max_degree << " (" << most_connected_actor << ")" << endl;
        }
        
        cout << "Memory optimization: Using integer IDs and deduplicated vectors" << endl;
        cout << "=============================\n" << endl;
    }
    
private:
    int countValidActors() const {
        int count = 0;
        for (const auto& actor : actors) {
            if (actor.id != -1) count++;
        }
        return count;
    }
    
    int countValidMovies() const {
        int count = 0;
        for (const auto& movie : movies) {
            if (movie.id != -1) count++;
        }
        return count;
    }
};

// Test with corrected small dataset
void testWithSmallDataset() {
    cout << "=== TESTING WITH SMALL DATASET ===" << endl;

    // name.basics.tsv
    ofstream name_file("test_names.tsv");
    name_file << "nconst\tprimaryName\tbirthYear\tdeathYear\tprimaryProfession\tknownForTitles\n";
    name_file << "nm1\tActor 1\t1950\t\\N\tactor\ttt1\n";
    name_file << "nm2\tActor 2\t1960\t\\N\tactor\ttt1,tt2\n";
    name_file << "nm3\tActor 3\t1970\t\\N\tactor\ttt2\n";
    name_file.close();

    // title.basics.tsv
    ofstream title_file("test_titles.tsv");
    title_file << "tconst\ttitleType\tprimaryTitle\toriginalTitle\tisAdult\tstartYear\tendYear\truntimeMinutes\tgenres\n";
    title_file << "tt1\tmovie\tMovie 1\tMovie 1\t0\t2000\t\\N\t120\tDrama\n";
    title_file << "tt2\tmovie\tMovie 2\tMovie 2\t0\t2010\t\\N\t110\tComedy\n";
    title_file.close();

    // title.principals.tsv
    ofstream principals_file("test_principals.tsv");
    principals_file << "tconst\tordering\tnconst\tcategory\tjob\tcharacters\n";
    principals_file << "tt1\t1\tnm1\tactor\t\\N\t[\"Char1\"]\n";
    principals_file << "tt1\t2\tnm2\tactor\t\\N\t[\"Char2\"]\n";
    principals_file << "tt2\t1\tnm2\tactor\t\\N\t[\"Char3\"]\n";
    principals_file << "tt2\t2\tnm3\tactor\t\\N\t[\"Char4\"]\n";
    principals_file.close();

    IMDBGraph test_graph;
    test_graph.loadActorsFromTSV("test_names.tsv");
    test_graph.loadMoviesFromTSV("test_titles.tsv");
    test_graph.loadPrincipalsFromTSV("test_principals.tsv");
    test_graph.buildCollaborationGraph();
    test_graph.printGraphStats();

    // Test shortest path
    int actor1 = test_graph.findActorByName("Actor 1");
    int actor3 = test_graph.findActorByName("Actor 3");
    if (actor1 != -1 && actor3 != -1) {
        auto path = test_graph.findShortestPath(actor1, actor3);
        cout << "Test path: ";
        for (int actor_id : path) cout << test_graph.getActorName(actor_id) << " -> ";
        cout << "END\n";
    }
}

// Main function
int main() {
    cout << "=== IMDB Social Network Analysis (Optimized) ===" << endl;
    
    // Test with small dataset first
    testWithSmallDataset();
    
    // Uncomment to run with full dataset
    /*
    IMDBGraph graph;
    
    cout << "Loading full IMDB dataset..." << endl;
    if (!graph.loadActorsFromTSV("name.basics.tsv")) {
        cerr << "Failed to load actors" << endl;
        return 1;
    }
    
    if (!graph.loadMoviesFromTSV("title.basics.tsv")) {
        cerr << "Failed to load movies" << endl;
        return 1;
    }
    
    if (!graph.loadPrincipalsFromTSV("title.principals.tsv")) {
        cerr << "Failed to load principals" << endl;
        return 1;
    }
    
    graph.buildCollaborationGraph();
    graph.printGraphStats();
    
    // Compute all centrality measures
    auto degree_centrality = graph.computeDegreeCentrality();
    auto betweenness_centrality = graph.computeBetweennessCentrality(500);
    auto closeness_centrality = graph.computeClosenessCentrality(500);
    
    // Get top actors
    auto top_degree = graph.getTopDegree(degree_centrality, 20);
    auto top_betweenness = graph.getTopCentrality(betweenness_centrality, 20);
    auto top_closeness = graph.getTopCentrality(closeness_centrality, 20);
    
    // Save results to CSV
    graph.saveDegreeResultsToCSV(top_degree, "degree_centrality.csv");
    graph.saveResultsToCSV(top_betweenness, "betweenness_centrality.csv", "Betweenness");
    graph.saveResultsToCSV(top_closeness, "closeness_centrality.csv", "Closeness");
    
    cout << "=== ANALYSIS COMPLETE ===" << endl;
    */
    
    return 0;
}