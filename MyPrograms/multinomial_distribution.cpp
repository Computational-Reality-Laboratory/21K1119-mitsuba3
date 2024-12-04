/*多項分布から乱数を生成*/
#include <iostream>
#include <random>


float get_rand() {
    // 乱数生成器（シードを指定可）
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0f, 1.0f);
    return dist(gen);
}

int main()
{
    // 試行回数
    int N = 100;
    // 確率ベクトル
    float p[4] = {1.0f/4.0f, 1.0f/4.0f, 1.0f/4.0f, 1.0f/4.0f};
    // 乱数生成時に使う乱数用
    float x = 0.0f;
    // 生成結果
    int result[4] = {0, 0, 0, 0};

    std::cout << "p = [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]" << std::endl;
    
    float P1 = p[0];
    float P2 = p[0] + p[1];
    float P3 = p[0] + p[1] + p[2];

    for(int i = 0; i < N; i++) {
        x = get_rand();
        if(x <= P1) {
            result[0] += 1;
        } 
        else if (x <= P2) {
            result[1] += 1;
        }
        else if (x <= P3) {
            result[2] += 1;
        }
        else {
            result[3] += 1;
        }
    }

    std::cout << "result = [" << result[0] << ", " << result[1] << ", " << result[2] << ", " << result[3] << "]" << std::endl;

}

