import java.util.Scanner;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Random;

public class QR_and_LU_decomposition {
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        Random random = new Random();
        int size = 1;

        System.out.println("Размер квадратной матрицы: (время выполнения LU), (время выполнения QR)");

        for (int k = 0; k < 4; k ++) {
            //System.out.print("Введите размерность квадратной матрицы: ");
            size *= 5;

            int size_row = size;
            int size_column = size;

            byte error;

            long duration1;
            long duration2;

            double[][] A = new double[size_row][size_column];
            double[][] Q = new double[size_row][size_column];
            double[][] R = new double[size_column][size_column];
            double[][] L = new double[size][size];
            double[][] U = new double[size][size];

            double [] b1 = new double[size];
            double [] b2 = new double[size];

            double [] x_with_LU_1 = new double[size];
            double [] x_with_LU_2 = new double[size];
            double [] x_with_QR_1 = new double[size];
            double [] x_with_QR_2 = new double[size];

            double [] y_with_LU_1 = new double[size];
            double [] y_with_LU_2 = new double[size];
            double [] y_with_QR_1 = new double[size];
            double [] y_with_QR_2 = new double[size];

            //System.out.println("Введите элементы матрицы построчно:");
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    //A[i][j] = scan.nextDouble();
                    A[i][j] = 15 * random.nextDouble();
                }
            }

            //System.out.println("Введите элементы вектора b1 через пробел:");
            for (int i = 0; i < size; i++) {
                //b1[i] = scan.nextDouble();
                b1[i] = 15 * random.nextDouble();
            }

            //System.out.println("Введите элементы вектора b2 через пробел:");
            for (int i = 0; i < size; i++) {
                //b2[i] = scan.nextDouble();
                b2[i] = 15 * random.nextDouble();
            }

            // Вызов LU-разложения после ввода данных
            long startTime = System.nanoTime();

            error = LU_decomposition(A, L, U, size);
            if (error == 0) {
                return ;
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        y_with_LU_1[i] = b1[i];
                    }
                    else {
                        y_with_LU_1[i] -= L[i][j-1] * y_with_LU_1[j-1];
                    }
                }
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        x_with_LU_1[(size - 1) - i] = y_with_LU_1[(size - 1) - i];
                    }
                    else {
                        x_with_LU_1[(size - 1) - i] -= U[(size - 1) - i][(size - 1) - (j - 1)] * x_with_LU_1[(size - 1) - (j - 1)];
                    }
                }
                x_with_LU_1[(size - 1) - i] /= U[(size - 1) - i][(size - 1) - i];
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        y_with_LU_2[i] = b2[i];
                    }
                    else {
                        y_with_LU_2[i] -= L[i][j-1] * y_with_LU_2[j-1];
                    }
                }
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        x_with_LU_2[(size - 1) - i] = y_with_LU_2[(size - 1) - i];
                    }
                    else {
                        x_with_LU_2[(size - 1) - i] -= U[(size - 1) - i][(size - 1) - (j - 1)] * x_with_LU_2[(size - 1) - (j - 1)];
                    }
                }
                x_with_LU_2[(size - 1) - i] /= U[(size - 1) - i][(size - 1) - i];
            }

            long endTime = System.nanoTime();
            duration1 = (endTime - startTime) / 1000;

            startTime = System.nanoTime();

            error = QR_decomposition(A, Q, R, size_row, size_column);
            if (error == 0) {
                return ;
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    y_with_QR_1[i] += Q[j][i] * b1[j];
                }
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        x_with_QR_1[(size - 1) - i] = y_with_QR_1[(size - 1) - i];
                    }
                    else {
                        x_with_QR_1[(size - 1) - i] -= R[(size - 1) - i][(size - 1) - (j - 1)] * x_with_QR_1[(size - 1) - (j - 1)];
                    }
                }
                x_with_QR_1[(size - 1) - i] /= R[(size - 1) - i][(size - 1) - i];
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    y_with_QR_2[i] += Q[j][i] * b2[j];
                }
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    if (j == 0) {
                        x_with_QR_2[(size - 1) - i] = y_with_QR_2[(size - 1) - i];
                    }
                    else {
                        x_with_QR_2[(size - 1) - i] -= R[(size - 1) - i][(size - 1) - (j - 1)] * x_with_QR_2[(size - 1) - (j - 1)];
                    }
                }
                x_with_QR_2[(size - 1) - i] /= R[(size - 1) - i][(size - 1) - i];
            }
            // Вывод результатов
            endTime = System.nanoTime();
            duration2 = (endTime - startTime) / 1000;

            // System.out.println();
            // System.out.println("Матрица U:");
            // Print_matrix(U, size, size);
            
            // System.out.println("Матрица L:");
            // Print_matrix(L, size, size);

            // System.out.println("Матрица Q:");
            // Print_matrix(Q, size_row, size_column);

            // System.out.println("Матрица R:");
            // Print_matrix(R, size_column, size_column);

            // Print_vector(x_with_LU_1, size);
            // Print_vector(x_with_LU_2, size);
            
            // Print_vector(x_with_QR_1, size);
            // Print_vector(x_with_QR_2, size);

            System.out.println(size +": "+ duration1 + " мкс, " + duration2 + " мкс");
        }

            scan.close();
    }

    public static byte LU_decomposition(double[][] A, double[][] L, double[][] U, int size) {

        // Копируем A в рабочую матрицу M
        double[][] M = new double[size][size];
        for (int i = 0; i < size; i++) {
            System.arraycopy(A[i], 0, M[i], 0, size);
        }

        for (int k = 0; k < size; k++) {
            if (M[k][k] == 0) {
                System.out.println("Данную матрицу разложить таким образом невозможно.");
                return 0;
            }

            // Формируем строку U[k][*]
            for (int j = k; j < size; j++) {
                U[k][j] = M[k][j];
            }

            // Формируем столбец L[*][k]
            for (int i = k; i < size; i++) {
                L[i][k] = M[i][k] / U[k][k];
            }

            // Обновляем рабочую матрицу M
            for (int i = k + 1; i < size; i++) {
                for (int j = k + 1; j < size; j++) {
                    M[i][j] = M[i][j] - L[i][k] * U[k][j];
                }
            }
        }

        return 1;
    }

    public static byte QR_decomposition(double[][] A, double[][] Q, double[][] R, int size_row, int size_column) {
        for (int i = 0; i < size_column; i++) {
            for (int j = 0; j <= i; j++) {
                if (j == 0) {
                    for (int k = 0; k < size_row; k++) {
                        Q[k][i] = A[k][i];
                    }
                }
                else {
                    for (int k = 0; k < size_row; k++) {
                        R[j-1][i] += Q[k][j-1] * A[k][i];
                    }
                    for (int k = 0; k < size_row; k++) {
                        Q[k][i] -= (R[j-1][i] * Q[k][j-1]);
                    }
                }
            }
            R[i][i] = Find_sqrt(Q, size_row, i);
            if (R[i][i] == 0) {
                System.out.println("Данная матрица не раскладывается на Q и R");
                return 0;
            }
            for (int k = 0; k < size_row; k++) {
                Q[k][i] /= R[i][i];
            }
        }

        return 1;
    }

    public static void Print_matrix(double[][] matrix, int size_row, int size_column) {
        for (int i = 0; i < size_row; i++) {
            for (int j = 0; j < size_column; j++) {
                System.out.printf("%.5f ", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();

    }

    public static void Print_vector(double[] vector, int size) {
        for (int i = 0; i < size; i++) {
            System.out.printf("%.5f ", vector[i]);
        }
        System.out.println();
    }

    // функция для подсчёта корня от суммы квадратов
    public static double Find_sqrt(double[][] matrix, int size_row, int column) {
        double sum = 0;
        for (int i = 0; i < size_row; i++) {
            BigDecimal bd = new BigDecimal(Double.toString(matrix[i][column] * matrix[i][column]));
            bd = bd.setScale(25, RoundingMode.HALF_UP); //Округление с точночтью до 25 знаков после запятой.
            sum += (bd.doubleValue());
        }
        return Math.sqrt(sum);
    }
}
