
template <typename T>
class hodlr 
{

private:
    int n, m;
    Matrix<T> a11, a22;
    LowRank<T> a21, a12;
};