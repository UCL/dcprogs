#include<type_traits>
int main()
{
  bool const a = std::is_move_constructible<int>::value;
  bool const a = std::is_move_assignable<int>::value;
  bool const a = std::is_trivially_move_assignable<int>::value;
  return 0;
}
