import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Random;
import java.util.Scanner;

public class QR_and_LU_decomposition {
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        Random random = new Random();
        int size = 1;

        System.out.println("Сравнение времени выполнения LU и QR разложения:");
        System.out.printf("%-25s %-25s %s%n", "Размер матрицы (n x n)", "LU Decomposition (мкс)", "QR Decomposition (мкс)");
        System.out.println("---------------------------------------------------------------------------------------------------");

        for (int k = 0; k < 4; k++) {
            size *= 5;

            int n = size;
            byte error;
            long duration_lu, duration_qr;

            double[][] A = new double[n][n];
            double[][] Q = new double[n][n];
            double[][] R = new double[n][n];
            double[][] L = new double[n][n];
            double[][] U = new double[n][n];

            double[] b1 = new double[n];
            double[] y_for_lu_1 = new double[n];
            double[] x_from_lu_1 = new double[n];

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = 15 * random.nextDouble();
                }
                b1[i] = 15 * random.nextDouble();
            }

            long startTime = System.nanoTime();

            error = LU_decomposition(A, L, U, n);
            if (error == 0) { 
                System.out.println("\nРазложение LU прервано.");
                break;
            }

            for (int i = 0; i < n; i++) {
                y_for_lu_1[i] = b1[i];
                for (int j = 0; j < i; j++) {
                    y_for_lu_1[i] -= L[i][j] * y_for_lu_1[j];
                }
            }

            for (int i = n - 1; i >= 0; i--) {
                x_from_lu_1[i] = y_for_lu_1[i];
                for (int j = i + 1; j < n; j++) {
                    x_from_lu_1[i] -= U[i][j] * x_from_lu_1[j];
                }
                x_from_lu_1[i] /= U[i][i];
            }

            long endTime = System.nanoTime();
            duration_lu = (endTime - startTime) / 1000;

            startTime = System.nanoTime();

            error = QR_decomposition(A, Q, R, n, n);
            if (error == 0) {
                System.out.println("\nРазложение QR прервано.");
                break;
            }

            double[] y_for_qr_1 = new double[n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    y_for_qr_1[i] += Q[j][i] * b1[j];
                }
            }

            double[] x_from_qr_1 = new double[n];
            for (int i = n - 1; i >= 0; i--) {
                x_from_qr_1[i] = y_for_qr_1[i];
                for (int j = i + 1; j < n; j++) {
                    x_from_qr_1[i] -= R[i][j] * x_from_qr_1[j];
                }
                x_from_qr_1[i] /= R[i][i];
            }

            endTime = System.nanoTime();
            duration_qr = (endTime - startTime) / 1000;

            System.out.printf("%-25d %-25d %d%n", n, duration_lu, duration_qr);
        }

        scan.close();
    }

    public static byte LU_decomposition(double[][] A, double[][] L, double[][] U, int size) {
        double[][] M = new double[size][size];
        for (int i = 0; i < size; i++) {
            System.arraycopy(A[i], 0, M[i], 0, size);
        }

        for (int k = 0; k < size; k++) {
            if (M[k][k] == 0) {
                System.out.println("Данную матрицу разложить таким образом (LU) невозможно.");
                return 0;
            }

            for (int j = k; j < size; j++) {
                U[k][j] = M[k][j];
            }

            for (int i = k; i < size; i++) {
                L[i][k] = M[i][k] / U[k][k];
            }

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
                } else {
                    for (int k = 0; k < size_row; k++) {
                        R[j-1][i] += Q[k][j-1] * A[k][i];
                    }
                    for (int k = 0; k < size_row; k++) {
                        Q[k][i] -= R[j-1][i] * Q[k][j-1];
                    }
                }
            }
            R[i][i] = Find_sqrt(Q, size_row, i);
            if (R[i][i] == 0) {
                System.out.println("Данная матрица не раскладывается на Q и R.");
                return 0;
            }
            for (int k = 0; k < size_row; k++) {
                Q[k][i] /= R[i][i];
            }
        }
        return 1;
    }

    public static double Find_sqrt(double[][] matrix, int size_row, int column) {
        double sum = 0;
        for (int i = 0; i < size_row; i++) {
            BigDecimal bd = new BigDecimal(Double.toString(matrix[i][column] * matrix[i][column]));
            bd = bd.setScale(25, RoundingMode.HALF_UP);
            sum += bd.doubleValue();
        }
        return Math.sqrt(sum);
    }
}
