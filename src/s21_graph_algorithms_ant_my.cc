#include "./graph/s21_graph.h"
#include "s21_graph_algorithms.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

const double ALPHA = 1.0;
const double BETA = 1.0;
const double RHO = 0.5;
const double QVAL = 100;
/*
ALPHA - коэффициент влиящий на определение следующей точки согласно количеству
феромонов BETA - коэффициент влиящий на определение следующей точки согласно
длины до точки RHO - коээфициент затухания феромонов QVAL - значение феромона
который откладывается
*/

// инциализация структуры муравьев
void GraphAlgorithms::initAnts(Variables &vars, std::vector<Ant> &antsVector) {
  int to, ant;
  /* Инициализация муравьев */
  to = 1;
  for (ant = 1; ant <= vars.Nvert; ant++) {
    // Распределяем муравьев по городам равномерно
    if (to > vars.Nvert)
      to = 1;
    antsVector[ant].curCity = to++;
    antsVector[ant].visited.assign(vars.Nvert + 1, false);
    antsVector[ant].pathIndex = 1;
    antsVector[ant].path.assign(vars.Nvert + 1, 1);
    antsVector[ant].path[0] = antsVector[ant].curCity;
    antsVector[ant].nextCity = 1;
    antsVector[ant].tourLength = 0.0;
    //  Помещаем исходный город, в котором находится муравей,
    // в список посещенных городов
    antsVector[ant].visited[antsVector[ant].curCity] = true;
  }
}

// заоплянет граф 0 по диагонали
void GraphAlgorithms::zeroGraph(Graph &graph, int N) {
  for (int i = 1; i <= N; i++) {
    graph.SetEdges(i, i, 0.0);
  }
}

// вычисление матрицы вероятностей для перехода в новую точку
void GraphAlgorithms::getProbability(Variables &vars, int N) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      if ((vars.grapPherom[i][j] > 0) || vars.graphDist.GetEdges(i, j) != 0) {
        vars.graphProbability[i][j] =
            (pow(vars.grapPherom[i][j], ALPHA) *
             pow((1.0 / vars.graphDist.GetEdges(i, j)), BETA));
      } else {
        vars.graphProbability[i][j] = 0;
      }
      // удалить возможно
      if (i == j)
        vars.graphProbability[i][j] = 0.0;
      else if (vars.graphProbability[i][j] < 0)
        vars.graphProbability[i][j] = 0.0;
    }
  }
}

// генератор случайных чисел от low до high
double GraphAlgorithms::randomValdouble(double low, double high) {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> distr(low, high);
  return distr(gen);
}

// выбор следующего города для муравья
int GraphAlgorithms::selectNextCity(int ant, int N, Variables &vars,
                                    std::vector<Ant> &antsVector) {
  int from, to;
  double denom = 0.0;
  double maxProb = 0;
  double temp = 0;
  bool flag = true;
  double counter = 0;

  /* Выбрать следующий город */
  from = antsVector[ant].curCity;
  /* Расчет знаменателя */
  for (to = 1; to <= N; to++) {
    if (antsVector[ant].visited[to] == 0) {
      temp = vars.graphProbability[from][to];
      denom += temp;
      if (temp > maxProb) {
        maxProb = temp;
      }
    }
  }
  maxProb = maxProb / denom;
  while (flag) {
    double p;
    to++;
    if (to > N)
      to = 1;
    if (antsVector[ant].visited[to] == 0) {
      p = vars.graphProbability[from][to] / denom;
      double zz = randomValdouble(0.0, maxProb);
      if (zz < p)
        flag = false;
    }
    if (counter > N * N) {
      if (vars.graphDist.GetEdges(from, to) > 0)
        flag = false;
    }
    counter++;
  }
  return to;
}

// запуск муравьев и определение их маршрута
int GraphAlgorithms::simulateAnts(std::vector<Ant> &antsVector, int N,
                                  Variables &vars) {
  int k;
  int moving = 0;
  for (k = 1; k <= N; k++) {
    /* Убедиться, что у муравья есть куда идти */
    if (antsVector[k].pathIndex < N) {
      antsVector[k].nextCity = selectNextCity(k, N, vars, antsVector);
      antsVector[k].visited[antsVector[k].nextCity] = 1;
      antsVector[k].path[antsVector[k].pathIndex++] = antsVector[k].nextCity;
      antsVector[k].tourLength += vars.graphDist.GetEdges(
          antsVector[k].curCity, antsVector[k].nextCity);

      /* Обработка окончания путешествия (из последнего города
       * в первый)
       */
      if (antsVector[k].pathIndex == N) {
        antsVector[k].tourLength += vars.graphDist.GetEdges(
            antsVector[k].path[N - 1], antsVector[k].path[0]);
      }
      antsVector[k].curCity = antsVector[k].nextCity;
      moving++;
    }
  }
  return moving;
}

// обновление матрицы ферментов
void GraphAlgorithms::updateTrails(std::vector<Ant> antsVector,
                                   std::vector<std::vector<double>> &grapPherom,
                                   int N, double init_pheromone) {
  int from, to, i, ant;
  // Испарение фермента
  for (from = 1; from <= N; from++) {
    for (to = 1; to <= N; to++) {
      if (from != to) {

        grapPherom[from][to] = (grapPherom[from][to] * (1.0 - RHO));
        if (grapPherom[from][to] <= 0)
          // grapPherom[from][to] = 0; // выбрать что удалить
          grapPherom[from][to] = init_pheromone; // выбрать что удалить
        // else
        //   grapPherom[from][to] = 0;
      }
    }
  }
  /* Нанесение нового фермента */
  /* Для пути каждого муравья */
  for (ant = 1; ant <= N; ant++) {
    /* Обновляем каждый шаг пути */
    for (i = 1; i <= N; i++) {
      if (i <= N - 1) {
        from = antsVector[ant].path[i];
        to = antsVector[ant].path[i + 1];
      } else {
        from = antsVector[ant].path[i];
        to = antsVector[ant].path[0];
      }
      grapPherom[from][to] =
          grapPherom[from][to] + QVAL / antsVector[ant].tourLength;
      grapPherom[to][from] = grapPherom[from][to];
    }
  }
  for (from = 1; from <= N; from++) {
    for (to = 1; to <= N; to++) {
      grapPherom[from][to] = grapPherom[from][to] * RHO;
    }
  }
}

// перезапуск муравьев
void GraphAlgorithms::restartAnts(std::vector<Ant> &antsVector, Variables &vars,
                                  int N, TsmResult &rez) {
  int ant, i, to = 1;
  for (ant = 1; ant <= N; ant++) {
    if (antsVector[ant].tourLength < vars.best) {
      vars.best = antsVector[ant].tourLength;
      vars.bestIndex = ant;
      std::vector<int> tmp(antsVector[ant].path.begin(),
                           antsVector[ant].path.end() - 1);
      rez.vertices = tmp;
    }
    antsVector[ant].nextCity = 1;
    antsVector[ant].tourLength = 0.0;
    for (i = 1; i <= N; i++) {
      antsVector[ant].visited[i] = 0;
      antsVector[ant].path[i] = 1;
    }
    if (to == N)
      to = 1;
    antsVector[ant].curCity = to++;
    antsVector[ant].pathIndex = 1;
    antsVector[ant].path[0] = antsVector[ant].curCity;
    antsVector[ant].visited[antsVector[ant].curCity] = 1;
  }
}

// решение задачи коммивояжера муравьиным алгоритмом
TsmResult GraphAlgorithms::SolveTravelingSalesmanProblem(Graph &graph) {
  TsmResult rez;
  Variables vars;
  vars.Nvert = (int)graph.GetVertices();
  vars.best = INF;
  vars.INIT_PHEROMONE = 1.0 / (double)vars.Nvert;
  vars.graphDist = graph;
  std::vector<std::vector<double>> temp(
      vars.Nvert + 1, std::vector<double>(vars.Nvert + 1, vars.INIT_PHEROMONE));
  vars.grapPherom = temp;
  std::vector<std::vector<double>> temp2(
      vars.Nvert + 1, std::vector<double>(vars.Nvert + 1, 0));
  vars.graphProbability = temp2;
  zeroGraph(vars.graphDist, vars.Nvert);
  std::vector<Ant> antsVector(vars.Nvert + 1);
  initAnts(vars, antsVector);
  getProbability(vars, vars.Nvert);
  int curTime = 0;
  double MAX_TIME = vars.Nvert * vars.Nvert * 25;
  while (curTime++ < MAX_TIME) {
    vars.curTime2 = curTime;
    if (simulateAnts(antsVector, vars.Nvert, vars) == 0) {
      updateTrails(antsVector, vars.grapPherom, vars.Nvert,
                   vars.INIT_PHEROMONE);
      getProbability(vars, vars.Nvert);
      if (curTime != MAX_TIME)
        restartAnts(antsVector, vars, vars.Nvert, rez);
    }
    // уменьшение феромонов для новых переходов если феромон обнулился
    vars.INIT_PHEROMONE *= 0.5;
  }
  return rez;
}
