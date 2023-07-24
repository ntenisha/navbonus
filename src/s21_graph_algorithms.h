#ifndef S21_GRAPH_ALGORITHMS_H
#define S21_GRAPH_ALGORITHMS_H

#include "./graph/s21_graph.h"
#include "./lib/list.h"
#include <vector>

const double INF = std::numeric_limits<int>::max();

struct TsmResult {
  // vector with the route you are looking for (with the vertex
  // traverse order).
  std::vector<int> vertices;
  // the length of this route
  double distance;
};

class NodeBB {
public:
  NodeBB() = default;
  std::vector<std::pair<int, int>> path;
  Graph graph;
  int cost;
  int vertex;
  int level;

private:
};

struct Ant {
  int curCity;
  int nextCity;
  std::vector<bool> visited;
  int pathIndex;
  std::vector<int> path;
  double tourLength;
};

struct Variables {
  Graph graphDist;
  std::vector<std::vector<double>> grapPherom;
  std::vector<std::vector<double>> graphProbability;
  int Nvert;
  double INIT_PHEROMONE;
  int curTime2;
};

class GraphAlgorithms {
public:
  std::vector<int> *DepthFirstSearch(const Graph &graph, int start_vertex);
  std::vector<int> BreadthFirstSearch(const Graph &graph, int start_vertex);
  double GetShortestPathBetweenVertices(const Graph &graph, int vertex1,
                                        int vertex2);
  // S21Matrix<double> *GetShortestPathsBetweenAllVertices(const Graph &graph);
  // S21Matrix<double> *GetLeastSpanningTree(const Graph &graph);
  // муравьи
  TsmResult SolveTravelingSalesmanProblem(Graph &graph);
  // границы и ветви
  TsmResult tspBoundnBounds(Graph &graph);
  // имитация отжига
  TsmResult annealingAlgorithm(Graph &graph);
  // эта функция участвует в консоли и расчитывает путь без кольца
  TsmResult getlength(std::vector<int> &path, Graph &graph);

  // todo  удалить определить это приватные или паблик методы
  bool isConnected(const Graph &graph);
  bool isDirected(const Graph &graph);

private:
  // границы и ветви
  void reduce_row(Graph &graph, std::vector<int> &row);
  void reduce_column(Graph &graph, std::vector<int> &col);
  int cost_calculation(Graph &graph);
  void infGraph(Graph &graph);
  // имитация отжига
  int randomVal(int low, int high);
  int randomValues(int low, int high);
  void getstablepath(TsmResult &stablerez, Graph graph);
  double ft_temp(double currentTemp, int Iteration);
  double distanceLenCycle(std::vector<int> &path, Graph &graph);
  // муравьи
  void initAnts(Variables &vars, std::vector<Ant> &antsVector);
  void zeroGraph(Graph &graph, int N);
  void getProbability(Variables &vars, int N);
  double randomValdouble(double low, double high);
  int selectNextCity(int ant, int N, Variables &vars,
                     std::vector<Ant> &antsVector);
  int simulateAnts(std::vector<Ant> &antsVector, int N, Variables &vars);
  void restartAnts(std::vector<Ant> &antsVector, int N, TsmResult &rez);

  void updateTrails(std::vector<Ant> antsVector,
                    std::vector<std::vector<double>> &grapPherom, int N,
                    double init_pheromone);
};

#endif // A2_SIMPLE_NAVIGATOR_SRC_S21_GRAPH_ALGORITHMS_H_
