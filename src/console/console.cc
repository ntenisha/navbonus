#include "console.h"

void ConsoleApp::PressAnyKey() {
  std::cout << "Press ENTER to continue . . . ";
  // std::cin.ignore();
}

bool ConsoleApp::getInput(int &num) {
  std::string input;
  std::cin >> input;
  ConsoleApp::ft_trim(input);

  try {
    num = std::stoi(input);
  } catch (const std::exception &e) {
    return false;
  }

  return true;
}

command ConsoleApp::getCommand() {
  while (true) {
    std::cout << std::endl
              << "\n-------------------------------------------------\n"
                 "Please, choose menu, then press key and  'enter': \n"
                 "#1: load graph from file \n"
                 "#2: breadth first search\n"
                 "#3: depth first search\n"
                 "#4: shortest path between a pair of vertices\n"
                 "#5: shortest path between all vertices\n"
                 "#6: least spanning tree\n"
                 "#7: solve travelling salesman problem\n"
                 "#8: compare time \n"
                 "#9: print dot\n"
                 "#10: print bnb \n"
                 "#11: print graph\n"
                 "#12: annealing simulation \n"
                 "#13: getlength \n"
                 "#0: quit \n"
                 "_______________________________________________________"
              << '\n';

    int number;
    if (!getInput(number)) {
      std::cout << "invalid input. Input should be of integer type. Try again."
                << std::endl;
      PressAnyKey();
      continue;
    }

    if (number < 0 || number > 13) {
      std::cout << number
                << " is out of bounds of expected command numbers. Try again"
                << std::endl;
      PressAnyKey();
      continue;
    }

    return (enum command)number;
  }
}

void ConsoleApp::ft_trim(std::string &str) {
  const std::string whiteSpaces = " \t\n\r\f\v";
  // Remove leading whitespace
  size_t first_non_space = str.find_first_not_of(whiteSpaces);
  str.erase(0, first_non_space);
  // Remove trailing whitespace
  size_t last_non_space = str.find_last_not_of(whiteSpaces);
  str.erase(last_non_space + 1);
  std::cout << str << std::endl;
}

void ConsoleApp::printVector(std::vector<int> vect1) {
  for (auto it : vect1) {
    std::cout << it << " ";
  }
  std::cout << "\n" << std::endl;
}

void ConsoleApp::printBFS(Graph &graph) {
  std::cout << "\nenter vertex number from 1 to " << graph.GetVertices()
            << " to start search " << std::endl;
  int vertex;
  GraphAlgorithms algor;
  std::cin >> vertex;
  try {
    std::vector<int> vecBFS;
    vecBFS = algor.BreadthFirstSearch(graph, vertex);
    printVector(vecBFS);
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::printDFS(Graph &graph) {
  std::cout << "\nenter vertex number from 1 to " << graph.GetVertices()
            << " to start search " << std::endl;
  int vertex;
  GraphAlgorithms algor;
  std::cin >> vertex;
  try {
    std::vector<int> *vecBFS;
    vecBFS = algor.DepthFirstSearch(graph, vertex);
    for (const int &vertex : *vecBFS) {
      std::cout << vertex << " ";
    }
    std::cout << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::printShortestBetweenTwo(Graph &graph) {
  std::cout << "\nenter vertex number from 1 to " << graph.GetVertices()
            << " to start search " << std::endl;
  int vertex1;
  std::cin >> vertex1;
  std::cout << "\nenter vertex number from 1 to " << graph.GetVertices()
            << " to start search " << std::endl;
  int vertex2;
  std::cin >> vertex2;
  GraphAlgorithms algor;
  try {
    double res;
    res = algor.GetShortestPathBetweenVertices(graph, vertex1, vertex2);
    std::cout << "\n" << res << "\n" << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

// добавить печать матрицы
void ConsoleApp::printShortestPathsBetweenAllVertices(Graph &graph) {
  GraphAlgorithms algor;

  try {
    S21Matrix<double> *res;
    res = algor.GetShortestPathsBetweenAllVertices(graph);

    // добавить печать матрицы
    // for (const double &vertex : *res) {
    //   std::cout << vertex << " ";
    // }
    // std::cout << std::endl;

  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

// добавить печать матрицы
void ConsoleApp::printLeastSpanningTree(Graph &graph) {
  GraphAlgorithms algor;

  try {
    S21Matrix<double> *res;
    res = algor.GetLeastSpanningTree(graph);

    // добавить печать матрицы
    // for (const double &vertex : *res) {
    //   std::cout << vertex << " ";
    // }
    // std::cout << std::endl;

  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::printBnB(Graph &graph) {
  GraphAlgorithms algor;
  try {
    std::vector<int> vecBFS;
    TsmResult Bnb = algor.tspBoundnBounds(graph);
    printVector(Bnb.vertices);
    std::cout << "distance  " << Bnb.distance << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::printAnnealing(Graph &graph) {
  GraphAlgorithms algor;
  try {
    std::vector<int> vecBFS;
    TsmResult Bnb = algor.annealingAlgorithm(graph);
    printVector(Bnb.vertices);
    std::cout << "distance  " << Bnb.distance << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::printTsmResult(Graph &graph) {
  GraphAlgorithms algor;
  try {
    TsmResult Bnb = algor.SolveTravelingSalesmanProblem(graph);
    printVector(Bnb.vertices);
    std::cout << "distance  " << Bnb.distance << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::algoGetLength(Graph &graph) {
  GraphAlgorithms algor;

  std::vector<int> vec;
  std::string input;
  int num = 0;

  // std::cout << "enter numbers" << std::endl;
  // std::cin >> input;
  // std::stringstream ss(input);
  // while (ss >> num) {
  //   vec.push_back(num);
  // }

  std::cout << "enter numbers" << std::endl;
  // TODO удалить дубль
  std::getline(std::cin, input);
  std::getline(std::cin, input);
  std::stringstream ss(input);
  while (ss >> num) {
    vec.push_back(num);
  }
  try {
    TsmResult Bnb = algor.getlength(vec, graph);
    printVector(Bnb.vertices);
    std::cout << "distance  " << Bnb.distance << std::endl;
  } catch (std::exception &e) {
    std::cerr << "\nError: " << e.what() << "\n";
  }
}

void ConsoleApp::executeCommand(command command, Graph &graph) {
  switch (command) {
  case loadGraph: {
    std::cout << "input path for the  filename to export " << std::endl;
    std::string filename;
    std::cin >> filename;

    graph.LoadGraphFromFile(filename);

    std::cout << "Graph has been successfully loaded!" << std::endl;
    break;
  }
  case breadthFirstSearch: {
    std::cout << "Breadth-first search" << std::endl;
    printBFS(graph);
    break;
  }
  case depthFirstSearch: {
    std::cout << "Depth-first search" << std::endl;
    printDFS(graph);
    break;
  }
  case shortestBetweenTwo: {
    std::cout << "GetShortestPathBetweenVertices" << std::endl;
    printShortestBetweenTwo(graph);
    break;
  }
  case shortestBetweenAll: {
    std::cout << "GetShortestPathsBetweenAllVertices" << std::endl;
    printShortestPathsBetweenAllVertices(graph);
    break;
  }
  case leastSpanningTree: {
    std::cout << "GetLeastSpanningTree" << std::endl;
    printLeastSpanningTree(graph);
    break;
  }
  case salesmanProblem: {
    std::cout << "SolveTravelingSalesmanProblem" << std::endl;
    printTsmResult(graph);
    break;
  }
  case tracetime: {
    std::cout << "bnb" << std::endl;
    printAntVsOtherAlgo(graph);
    break;
  }
  case printDot: {
    std::cout << "print Dot file" << std::endl;
    break;
  }
  case quit: {
    std::cout << "Exit / Quit" << std::endl;
    exit(0);
    break;
  }
  case boundnbounds: {
    std::cout << "bnb" << std::endl;
    printBnB(graph);
    break;
  }
  case printGraph1: {
    std::cout << "print graph" << std::endl;
    graph.PrintAdjacencyMatrix();
    break;
  }

  case annealing: {
    std::cout << "annealing" << std::endl;
    printAnnealing(graph);
    break;
  }
  case getLength: {
    std::cout << "getLength" << std::endl;
    algoGetLength(graph);
    break;
  }
  }
}

void ConsoleApp::start() {
  /*	добавить

      инициализация графа,матрицы и прочего

  */
  Graph graph;

  while (true) {
    command command = getCommand();
    // if (command == quit) {
    //   return;
    // }
    try {
      /*	заменить

            executeCommand(command, graph);

      */
      executeCommand(command, graph);
      PressAnyKey();
    } catch (std::invalid_argument &ex) {
      std::cout << "invalid argument: " << ex.what() << std::endl;
      PressAnyKey();
    } catch (std::exception &ex) {
      std::cout << "unexpected exception: " << ex.what() << std::endl;
      PressAnyKey();
    }
  }
}

void ConsoleApp::printAntVsOtherAlgo(Graph &graph) {
  try {
    std::cout << "\nenter how many times to keep track of the time (1 - 1000)"
              << std::endl;
    int N;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      GraphAlgorithms algor;
      double timeAnt = 0.0;
      double timeAnneAlingAlgorithm = 0.0;
      double timeBranchAndBound = 0.0;
      TsmResult antRez;
      TsmResult anneAlingRez;
      TsmResult branchAndBoundRez;

      clock_t start, end;

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          antRez = algor.SolveTravelingSalesmanProblem(graph);
        } else {
          TsmResult tmp = algor.SolveTravelingSalesmanProblem(graph);
          if (tmp.distance < antRez.distance) {
            antRez = tmp;
          }
        }
      }
      end = clock();
      timeAnt = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          anneAlingRez = algor.annealingAlgorithm(graph);
        } else {
          TsmResult tmp = algor.annealingAlgorithm(graph);
          if (tmp.distance < anneAlingRez.distance) {
            anneAlingRez = tmp;
          }
        }
      }
      end = clock();
      timeAnneAlingAlgorithm = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      start = clock();
      for (int i = 0; i < N; ++i) {
        if (i == 0) {
          branchAndBoundRez = algor.tspBoundnBounds(graph);
        } else {
          TsmResult tmp = algor.tspBoundnBounds(graph);
          if (tmp.distance < branchAndBoundRez.distance) {
            branchAndBoundRez = tmp;
          }
        }
      }
      end = clock();
      timeBranchAndBound = ((double)end - start) / ((double)CLOCKS_PER_SEC);

      std::cout << "1) Ant Algoritm\n";
      std::cout << "REPS: " << N << std::endl;
      std::cout << "TIME: " << timeAnt << "\nDISTANCE: " << antRez.distance
                << std::endl;

      std::cout << "2) Annealing Algoritm\n";
      std::cout << "REPS: " << N << std::endl;
      std::cout << "TIME: " << timeAnneAlingAlgorithm
                << "\nDISTANCE: " << anneAlingRez.distance << std::endl;

      std::cout << "3) Branch and Bound Algoritm\n";
      std::cout << "REPS: " << N << std::endl;
      std::cout << "TIME: " << timeBranchAndBound
                << "\nDISTANCE: " << branchAndBoundRez.distance << std::endl;

    } else {
      std::cout << "\nINCORRECT COUNT OF REPS!" << std::endl;
    }
  } catch (std::exception &exceptionText) {
    std::cerr << "\nError: " << exceptionText.what() << "\n";
  }
}
