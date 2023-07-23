#include "./graph/s21_graph.h"
#include "s21_graph_algorithms.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <stack>
#include <vector>

// создание нода
NodeBB *newNode(Graph &graph, std::vector<std::pair<int, int>> const &path,
                int level, int i, int j) {
  NodeBB *node = new NodeBB();
  node->path = path;
  if (level != 1)
    node->path.push_back(std::make_pair(i, j));
  node->graph = graph; // копируем граф

  for (size_t k = 1; level != 1 && k <= graph.GetVertices(); k++) {
    node->graph.SetEdges(i, k, INF);
    node->graph.SetEdges(k, j, INF);
  }

  node->graph.SetEdges(j, 1, INF);
  node->level = level;
  node->vertex = j;
  return node;
}

// редуцирование строк матрицы
void GraphAlgorithms::reduce_row(Graph &graph, std::vector<int> &row) {
  int N = row.size();
  std::fill_n(row.begin(), N, INF);

  // находим минимальное значение в каждой строке
  for (int i = 1; i < N; i++)
    for (int j = 1; j < N; j++)
      if (graph.GetEdges(i, j) < row[i])
        row[i] = graph.GetEdges(i, j);

  // вычитаем минимальное значение из каждого элемента строки
  for (int i = 1; i < N; i++)
    for (int j = 1; j < N; j++)
      if (graph.GetEdges(i, j) != INF && row[i] != INF)
        graph.SetEdges(i, j, graph.GetEdges(i, j) - row[i]);
}

// редуцирование столбцов матрицы
void GraphAlgorithms::reduce_column(Graph &graph, std::vector<int> &col) {
  int N = col.size();
  std::fill_n(col.begin(), N, INF);

  // находим минимальное значение в каждом столбце
  for (int i = 1; i < N; i++)
    for (int j = 1; j < N; j++)
      if (graph.GetEdges(i, j) < col[j])
        col[j] = graph.GetEdges(i, j);

  // вычитаем минимальное значение из каждого элемента столбца
  for (int i = 1; i < N; i++)
    for (int j = 1; j < N; j++)
      if (graph.GetEdges(i, j) != INF && col[j] != INF)
        graph.SetEdges(i, j, graph.GetEdges(i, j) - col[j]);
}

// вычисление стоимости пути через редуцированную матрицу
int GraphAlgorithms::cost_calculation(Graph &graph) {
  int cost = 0;
  int N = graph.GetVertices() + 1;

  std::vector<int> row(N);
  reduce_row(graph, row);
  std::vector<int> col(N);
  reduce_column(graph, col);

  // суммируем значения минимальных элементов строк и столбцов
  for (int i = 1; i < N; i++) {
    cost += (row[i] != INF) ? row[i] : 0, cost += (col[i] != INF) ? col[i] : 0;
  }
  return cost;
}

// функционатор для очереди с приоритеттом
class comp {
public:
  bool operator()(const NodeBB *lhs, const NodeBB *rhs) const {
    return lhs->cost > rhs->cost;
  }
};

// заоплянет граф inf по диагонали
void GraphAlgorithms::infGraph(Graph &graph) {
  int N = graph.GetVertices() + 1;
  // удалить
  //  int i = 1;
  //  while (i < N) {
  //    graph.SetEdges(i, i, INF);
  //    i++;
  //  }
  for (int i = 1; i < N; i++) {
    graph.SetEdges(i, i, INF);
  }
}

// решение задачи коммивояжера методом ветвей и границ
TsmResult GraphAlgorithms::tspBoundnBounds(Graph &graph) {
  std::priority_queue<NodeBB *, std::vector<NodeBB *>, comp> pq;
  std::vector<std::pair<int, int>> v;
  //  Graph graphNew = graph; // если все норм то удалить, если нет то заменить
  //  на нижнюю
  Graph graphNew(graph);
  int N = graph.GetVertices() + 1;
  infGraph(graphNew);
  NodeBB *root = newNode(graphNew, v, 1, 1, 1);
  root->cost = cost_calculation(root->graph);
  pq.push(root);
  TsmResult result;
  // добавляем переменную для хранения лучшей найденной стоимости
  double best_cost = INF;
  while (!pq.empty()) {
    NodeBB *min = pq.top();
    pq.pop();
    int i = min->vertex;
    if (min->level == N - 1) {
      min->path.push_back(std::make_pair(i, 0));
      result.vertices.resize(N - 1);
      for (size_t k = 0; k < min->path.size(); k++) {
        result.vertices[k] = min->path[k].first;
      }
      result.distance = min->cost;
      // обновляем лучшую найденную стоимость
      best_cost = min->cost;
      delete min;
      return result;
    }

    // проверяем превышение текущей лучшей найденной стоимости
    if (min->cost >= best_cost) {
      // если превышение, удаляем узел и возвращаемся к предыдущему узлу
      delete min;
      continue;
    }
    for (int j = 1; j < N; j++) {
      if (min->graph.GetEdges(i, j) != INF) {
        NodeBB *child = newNode(min->graph, min->path, min->level + 1, i, j);
        child->cost = min->cost + min->graph.GetEdges(i, j) +
                      cost_calculation(child->graph);
        // добавляем только те узлы, которые могут привести к лучшему решению
        if (child->cost < best_cost)
          pq.push(child);
        // если превышение, удаляем узел
        else
          delete child;
      }
    }
    delete min;
  }
  return result;
}

// вычисляет путь с замыканием последней точки с 1
double GraphAlgorithms::distanceLenCycle(std::vector<int> &path, Graph &graph) {
  double dis = 0;
  for (size_t i = 1; i < path.size(); i++) {
    if (graph.GetEdges(path[i - 1], path[i]) > 0 &&
        graph.GetEdges(path[i - 1], path[i]) != INF) {
      dis += graph.GetEdges(path[i - 1], path[i]);
    } else {
      dis = INF;
      break;
    }
  }
  if (graph.GetEdges(path[path.size() - 1], path[0]) > 0 &&
      graph.GetEdges(path[path.size() - 1], path[0]) != INF && dis != INF) {
    dis += graph.GetEdges(path[path.size() - 1], path[0]);
  }
  return dis;
}

// вычисляет путь по точкам
TsmResult GraphAlgorithms::getlength(std::vector<int> &path, Graph &graph) {
  TsmResult res{{}, 0};
  res.vertices.push_back(0);
  for (size_t i = 1; i < path.size(); i++) {
    res.distance += graph.GetEdges(path[i - 1], path[i]);
    res.vertices.push_back(graph.GetEdges(path[i - 1], path[i]));
  }
  return res;
}

// в одно из вариантов не запускается на винде
// рандомизатор для int
int GraphAlgorithms::randomVal(int low, int high) {
  std::random_device rd;
  std::default_random_engine gen(rd());
  //  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> distr(low, high);
  return distr(gen);
}

// инциализация 1 пути путем рандомного заполнения
std::vector<int> randomVector(int low, int hgh) {
  std::vector<int> numbers(hgh - low + 1);
  // заполнение вектора числами в заданном диапазоне
  for (int i = low; i <= hgh; i++) {
    numbers[i - low] = i;
  }
  // перемешивание вектора в случайном порядке
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(numbers.begin(), numbers.end(), g);

  return numbers;
}

// рандомное перставление в пути 2 точек
void swapVector(std::vector<int> &myVector) {
  // генерируем два случайных индекса элементов в векторе
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, myVector.size() - 1);
  int index1 = dis(gen);
  int index2 = dis(gen);
  // меняем местами элементы с указанными индексами
  std::swap(myVector[index1], myVector[index2]);
}

// вычисление пути пока не найден будет путь ,на случай если добавиться INF в
// путь
void GraphAlgorithms::getstablepath(TsmResult &stablerez, Graph graph) {
  stablerez.distance = distanceLenCycle(stablerez.vertices, graph);
  while (stablerez.distance == INF) {
    stablerez.distance = distanceLenCycle(stablerez.vertices, graph);
    if (stablerez.distance == INF) {
      swapVector(stablerez.vertices);
    }
  }
}

// функция уменьшения температуры
double GraphAlgorithms::ft_temp(double currentTemp, int Iteration) {
  double newTemp;
  // закон больцмана
  // newTemp = currentTemp / (std::log(1 + Iteration));
  // закон коши
  // newTemp = currentTemp / (1 + Iteration);
  // верхние 2 слишком быстро уменьшают поэтому взял такую
  newTemp = currentTemp * Iteration / (1 + Iteration);
  return newTemp;
}

// решение задачи комивояжера путем имитации отжига
TsmResult GraphAlgorithms::annealingAlgorithm(Graph &graph) {
  TsmResult tmprez{{}, INF};
  TsmResult rez{{}, INF};
  // ввел minrez чтобы наверняка был самый короткий из найденых путей , иначе
  // будет будет больше путь
  TsmResult minrez{{}, INF};
  Graph graphNew = graph;
  double currentTemp = 500;
  double endTemp = 0.009;
  int Iteration = 0;
  int N = (int)graphNew.GetVertices();
  // ввел переменную maxIteration потому что по моей функции температура слишком
  // долго уменьшается
  double maxIteration = N * N * 25;
  infGraph(graphNew);
  rez.vertices = randomVector(1, N);
  getstablepath(rez, graphNew);

  tmprez = rez;
  minrez = rez;
  while (currentTemp > endTemp && Iteration < maxIteration) {
    swapVector(tmprez.vertices);
    getstablepath(tmprez, graphNew);
    if (tmprez.distance < rez.distance) {
      rez = tmprez;
      if (rez.distance < minrez.distance)
        minrez = rez;
    } else if ((int)(exp((rez.distance - tmprez.distance) / currentTemp) *
                     100) > (randomVal(0, 100))) {
      rez = tmprez;
    }
    Iteration++;
    currentTemp = ft_temp(currentTemp, Iteration);
  }
  return minrez;
}
