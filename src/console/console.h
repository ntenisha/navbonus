#ifndef CONSOLE_H
#define CONSOLE_H

#include "../graph/s21_graph.h"
#include "../s21_graph_algorithms.h"
#include <iostream>
#include <string>
#include <vector>

enum command {
  quit = 0,
  loadGraph = 1,
  breadthFirstSearch = 2,
  depthFirstSearch = 3,
  shortestBetweenTwo = 4,
  shortestBetweenAll = 5,
  leastSpanningTree = 6,
  salesmanProblem = 7,
  tracetime = 8,
  printDot = 9,
  boundnbounds = 10,
  printGraph1 = 11,
  annealing = 12,
  getLength = 13
};

class ConsoleApp {
public:
  void start();

private:
  void PressAnyKey();
  bool getInput(int &num);
  command getCommand();
  void executeCommand(command command, Graph &graph);
  void ft_trim(std::string &s);
  void printVector(std::vector<int> vect1);


  void printBFS(Graph &graph);
  void printDFS(Graph &graph);
  void printShortestBetweenTwo(Graph &graph);
  void printShortestPathsBetweenAllVertices(Graph &graph);
  void printLeastSpanningTree(Graph &graph);
  void printTsmResult(Graph &graph);
  void printAntVsOtherAlgo(Graph &graph);

  void printBnB(Graph &graph);
  void printGraph(Graph &graph);
  void printAnnealing(Graph &graph);
  void algoGetLength(Graph &graph);
};

#endif
