#if __clang_major__ < 3 || (__clang_major__ == 3 && __clang_minor__ < 1)
        #error Minimum Clang version supported is 3.1.
#endif

int main()
{
    return 0;
}
