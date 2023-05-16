#include <bits/stdc++.h>

#include <optional>

using namespace std;

const int kSize = 16;
const int kInf = 1000000;

class Vertex {
 public:
  enum class Directions { Up, Down, Left, Right };
  size_t matr[16];
  explicit Vertex(size_t a) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        matr[i * 4 + j] = a & 15ull;
        a >>= 4;
      }
    }
  }

  int FindZero() {
    int x;
    for (int i = 0; i < kSize; ++i) {
      if (matr[i] == 15) {
        x = i;
        break;
      }
    }
    return x;
  }

  std::optional<Vertex> GoUp() {
    int index = FindZero();
    if (index / 4 == 3) {
      return std::nullopt;
    }
    Vertex temp = *this;
    swap(temp.matr[index], temp.matr[index + 4]);
    return {temp};
  }

  std::optional<Vertex> GoDown() {
    int index = FindZero();
    if (index / 4 == 0) {
      return std::nullopt;
    }
    Vertex temp = *this;
    swap(temp.matr[index], temp.matr[index - 4]);
    return {temp};
  }

  std::optional<Vertex> GoLeft() {
    int index = FindZero();
    if (index % 4 == 3) {
      return std::nullopt;
    }
    Vertex temp = *this;
    swap(temp.matr[index], temp.matr[index + 1]);
    return {temp};
  }

  std::optional<Vertex> GoRight() {
    int index = FindZero();
    if (index % 4 == 0) {
      return std::nullopt;
    }
    Vertex temp = *this;
    swap(temp.matr[index], temp.matr[index - 1]);
    return {temp};
  }

  size_t LeftHeuristic() {
    size_t res = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (matr[i * 4 + j] == 15) {
          continue;
        }
        int usual_x = (matr[i * 4 + j]) / 4;
        int usual_y = (matr[i * 4 + j]) % 4;
        res += abs(i - usual_x) + abs(j - usual_y);
      }
    }
    for (int line = 0; line < 4; ++line) {
      for (int l = 0; l < 3; ++l) {
        for (int r = l + 1; r < 4; ++r) {
          if (matr[line * 4 + l] / 4 == line && matr[line * 4 + r] / 4 == line &&
              matr[line * 4 + l] > matr[line * 4 + r] && ((matr[line * 4 + l] | matr[line * 4 + r]) ^ 15)) {
            res += 2;
          }
        }
      }
    }
    for (int row = 0; row < 4; ++row) {
      for (int l = 0; l < 3; ++l) {
        for (int r = l + 1; r < 4; ++r) {
          if (matr[l * 4 + row] % 4 == row && matr[r + row % 4] % 4 == row &&
              matr[l * 4 + row] > matr[r * 4 + row] && ((matr[l * 4 + row] | matr[r * 4 + row]) ^ 15)) {
            res += 2;
          }
        }
      }
    }
    return (res * 111) / 100;
  }

  size_t RightHeurestic(Vertex& left) {
    int pos[16];
    for (int i = 0; i < 16; ++i) {
      pos[left.matr[i]] = i;
    }

    size_t res = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (matr[i * 4 + j] == 15) {
          continue;
        }
        int usual_x = (pos[matr[i * 4 + j]]) / 4;
        int usual_y = (pos[matr[i * 4 + j]]) % 4;
        res += abs(i - usual_x) + abs(j - usual_y);
      }
    }
    for (int line = 0; line < 4; ++line) {
      for (int l = 0; l < 3; ++l) {
        for (int r = l + 1; r < 4; ++r) {
          if (pos[matr[line * 4 + l]] / 4 == line && pos[matr[line * 4 + r]] / 4 == line &&
              pos[matr[line * 4 + l]] % 4 > pos[matr[line * 4 + r]] % 4 &&
              ((matr[line * 4 + l] | matr[line * 4 + r]) ^ 15)) {
            res += 2;
          }
        }
      }
    }
    for (int row = 0; row < 4; ++row) {
      for (int l = 0; l < 3; ++l) {
        for (int r = l + 1; r < 4; ++r) {
          if (pos[matr[l * 4 + row]] % 4 == row && pos[matr[r + row % 4]] % 4 == row &&
              pos[matr[l * 4 + row]] / 4 > pos[matr[r * 4 + row]] / 4 &&
              ((matr[l * 4 + row] | matr[r * 4 + row]) ^ 15)) {
            res += 2;
          }
        }
      }
    }
    return (res * 111) / 100;
  }

  vector<pair<size_t, Directions>> FindNext() {
    vector<pair<size_t, Directions>> ans;
    std::optional<Vertex> up_ans = GoUp();
    if (up_ans != nullopt) {
      ans.push_back(make_pair(up_ans.value().to_int(), Directions::Up));
    }

    std::optional<Vertex> down_ans = GoDown();
    if (down_ans != nullopt) {
      ans.push_back(make_pair(down_ans.value().to_int(), Directions::Down));
    }

    std::optional<Vertex> left_ans = GoLeft();
    if (left_ans != nullopt) {
      ans.push_back(make_pair(left_ans.value().to_int(), Directions::Left));
    }

    std::optional<Vertex> right_ans = GoRight();
    if (right_ans != nullopt) {
      ans.push_back(make_pair(right_ans.value().to_int(), Directions::Right));
    }
    return ans;
  }

  size_t to_int() {
    size_t res = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        res |= matr[i * 4 + j] * 1ull << (i * 4 + j) * 4;
      }
    }
    return res;
  }
};

bool Check(vector<size_t>& start) {
  int inv = 0;
  for (int i = 0; i < 16; ++i)
    if (start[i] != 15)
      for (int j = 0; j < i; ++j)
        if (start[j] > start[i] && start[j] != 15) ++inv;
  for (int i = 0; i < 16; ++i)
    if (start[i] == 15) inv += 1 + i / 4;

  return inv % 2 == 0;
}

size_t AStar(size_t start, size_t end, unordered_map<size_t, size_t>& left_parents,
             unordered_map<size_t, size_t>& right_parents) {
  priority_queue<pair<size_t, size_t>, std::vector<pair<size_t, size_t>>, greater<>> left_q;
  priority_queue<pair<size_t, size_t>, std::vector<pair<size_t, size_t>>, greater<>> right_q;
  unordered_map<size_t, size_t> g_left;
  unordered_map<size_t, size_t> g_right;
  unordered_set<size_t> left_used, right_used;
  g_left[start] = 0;
  g_right[end] = 0;
  size_t cur_min_dist = kInf;
  left_q.push({Vertex(start).LeftHeuristic(), start});
  size_t cur_mid;
  Vertex vertex_start(start);
  right_q.push({Vertex(end).RightHeurestic(vertex_start), end});
  left_used.insert(start);
  right_used.insert(end);
  while (cur_min_dist > max(left_q.top().first, right_q.top().first)) {
    size_t tmp = left_q.top().second;
    left_q.pop();
    if (right_used.find(tmp) != right_used.end() && g_left[tmp] + g_right[tmp] < cur_min_dist) {
      cur_min_dist = g_left[tmp] + g_right[tmp];
      cur_mid = tmp;
    }
    left_used.insert(tmp);
    vector<pair<size_t, Vertex::Directions>> next_direcitons = Vertex(tmp).FindNext();
    size_t dist = g_left[tmp];
    for (auto v : next_direcitons) {
      size_t next_weight = dist + 1;
      if (g_left.find(v.first) == g_left.end() || next_weight < g_left[v.first]) {
        left_parents[v.first] = tmp;
        g_left[v.first] = next_weight;
        left_q.push({next_weight + Vertex(v.first).LeftHeuristic(), v.first});
      }
    }

    tmp = right_q.top().second;
    right_q.pop();
    if (left_used.find(tmp) != left_used.end() && g_left[tmp] + g_right[tmp] < cur_min_dist) {
      cur_min_dist = g_left[tmp] + g_right[tmp];
      cur_mid = tmp;
    }
    right_used.insert(tmp);
    next_direcitons = Vertex(tmp).FindNext();
    dist = g_right[tmp];
    for (auto v : next_direcitons) {
      size_t next_weight = dist + 1;
      if (g_right.find(v.first) == g_right.end() || next_weight < g_right[v.first]) {
        right_parents[v.first] = tmp;
        g_right[v.first] = next_weight;
        right_q.push({next_weight + Vertex(v.first).RightHeurestic(vertex_start), v.first});
      }
    }
  }
  return cur_mid;
}

string PrintAns(size_t start, size_t end, size_t mid, unordered_map<size_t, size_t>& left_parents,
                unordered_map<size_t, size_t>& right_parents) {
  vector<Vertex::Directions> ans;
  size_t tmp = mid;
  string ans_string;
  while (tmp != start) {
    size_t prev = left_parents[tmp];
    vector<pair<size_t, Vertex::Directions>> next_direcitons = Vertex(tmp).FindNext();
    for (auto v : next_direcitons) {
      if (v.first == prev) {
        ans.push_back(v.second);
        break;
      }
    }
    tmp = prev;
  }
  reverse(ans.begin(), ans.end());
  for (size_t i = 0; i < ans.size(); ++i) {
    if (ans[i] == Vertex::Directions::Up) {
      ans_string += 'D';
    } else if (ans[i] == Vertex::Directions::Down) {
      ans_string += 'U';
    } else if (ans[i] == Vertex::Directions::Left) {
      ans_string += 'R';
    } else if (ans[i] == Vertex::Directions::Right) {
      ans_string += 'L';
    }
  }
  ans = {};
  tmp = mid;
  while (tmp != end) {
    size_t prev = right_parents[tmp];
    vector<pair<size_t, Vertex::Directions>> next_direcitons = Vertex(tmp).FindNext();
    for (auto v : next_direcitons) {
      if (v.first == prev) {
        ans.push_back(v.second);
        break;
      }
    }
    tmp = prev;
  }
  for (size_t i = 0; i < ans.size(); ++i) {
    if (ans[i] == Vertex::Directions::Up) {
      ans_string += 'U';
    } else if (ans[i] == Vertex::Directions::Down) {
      ans_string += 'D';
    } else if (ans[i] == Vertex::Directions::Left) {
      ans_string += 'L';
    } else if (ans[i] == Vertex::Directions::Right) {
      ans_string += 'R';
    }
  }
  return ans_string;
}

signed main() {
  cin.tie(0);
  ios_base::sync_with_stdio(false);
  vector<size_t> matrix(kSize);
  size_t start = 0;
  for (size_t i = 0; i < kSize; ++i) {
    cin >> matrix[i];
    if (matrix[i] == 0)
      matrix[i] = 15;
    else
      --matrix[i];
    start |= matrix[i] << i * 4;
  }
  if (!Check(matrix)) {
    cout << -1;
    return 0;
  }
  size_t end = 0xfedcba9876543210ull;
  unordered_map<size_t, size_t> left_parents, right_parents;
  left_parents.reserve(1000000);
  right_parents.reserve(1000000);
  size_t tmp = AStar(start, end, left_parents, right_parents);
  string ans = PrintAns(start, end, tmp, left_parents, right_parents);
  cout << ans.size() << '\n' << ans;
  return 0;
}
