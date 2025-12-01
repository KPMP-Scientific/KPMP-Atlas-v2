using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Common_functions
{
    class Correlation_coefficient_class
    {
        public double Get_pearson_correlation_coefficient_numerically_instable(double[] array_x, double[] array_y)
        {
            int array_x_length = array_x.Length;
            int array_y_length = array_y.Length;
            if (array_x_length != array_y_length)
            {
                throw new Exception();
            }

            double[] array_xy = new double[array_x_length];
            double[] array_xp2 = new double[array_x_length];
            double[] array_yp2 = new double[array_x_length];
            double sum_x = 0;
            double sum_y = 0;
            double sum_xy = 0;
            double sum_xpow2 = 0;
            double sum_ypow2 = 0;
            for (int indexXY = 0; indexXY < array_x_length; indexXY++)
            {
                array_xy[indexXY] = array_x[indexXY] * array_y[indexXY];
            }
            for (int indexXP2 = 0; indexXP2 < array_x_length; indexXP2++)
            {
                array_xp2[indexXP2] = Math.Pow(array_x[indexXP2], 2.0);
            }
            for (int indexYP2 = 0; indexYP2 < array_x_length; indexYP2++)
            {
                array_yp2[indexYP2] = Math.Pow(array_y[indexYP2], 2.0);
            }
            for (int indexYP2 = 0; indexYP2 < array_y_length; indexYP2++)
            {
                array_yp2[indexYP2] = Math.Pow(array_y[indexYP2], 2.0);
            }
            for (int index_x = 0; index_x < array_x_length; index_x++)
            {
                sum_x += array_x[index_x];
            }
            for (int index_y = 0; index_y < array_y_length; index_y++)
            {
                sum_y += array_y[index_y];
            }
            for (int index_xy = 0; index_xy < array_x_length; index_xy++)
            {
                sum_xy += array_xy[index_xy];
            }
            for (int indexXpow2 = 0; indexXpow2 < array_x_length; indexXpow2++)
            {
                sum_xpow2 += array_xp2[indexXpow2];
            }
            for (int indexYpow2 = 0; indexYpow2 < array_y_length; indexYpow2++)
            {
                sum_ypow2 += array_yp2[indexYpow2];
            }
            double Ex2 = Math.Pow(sum_x, 2.00);
            double Ey2 = Math.Pow(sum_y, 2.00);
            double Correl = (array_x_length * sum_xy - sum_x * sum_y) / Math.Sqrt((array_x_length * sum_xpow2 - Ex2) * (array_x_length * sum_ypow2 - Ey2));
            return Correl;
        }

        public double Get_pearson_correlation_coefficient_stable(double[] array_x, double[] array_y)
        {
            int n = array_x.Length;
            if (n != array_y.Length) { throw new ArgumentException("Length mismatch"); }

            double mean_x = array_x.Average();
            double mean_y = array_y.Average();

            double sum_xy = 0;
            double sum_x2 = 0;
            double sum_y2 = 0;

            //double max_double = 1E-10 * Double.MaxValue;
            //double quotient_xy = 1;
            //double quotient_x2 = 1;
            //double quotient_y2 = 1;


            for (int i = 0; i < n; i++)
            {
                double dx = array_x[i] - mean_x;
                double dy = array_y[i] - mean_y;
                //if ((sum_xy >= max_double)&&(quotient_xy==1))
                //{
                //quotient_xy = max_double;
                //sum_xy /= quotient_xy;
                //}
                //if ((sum_x2 >= max_double)&& (quotient_x2 == 1))
                //{
                //quotient_x2 = max_double;
                //sum_x2 /= quotient_x2;
                //}
                //if ((sum_y2 >= max_double)&& (quotient_y2 == 1))
                //{
                //quotient_y2 = max_double;
                //sum_y2 /= quotient_y2;
                //}

                sum_xy += dx * dy;// quotient_xy;
                sum_x2 += dx * dx;// quotient_x2;
                sum_y2 += dy * dy;// quotient_y2;
            }

            double denominator = Math.Sqrt(sum_x2 * sum_y2);
            if (Double.IsInfinity(denominator))
            {
                denominator = Math.Sqrt(sum_x2) * Math.Sqrt(sum_y2);
            }

            if (denominator == 0) { throw new Exception("Boundary problem in Get_pearson_correlation_coefficient_stable"); }
            return (sum_xy / denominator);
        }

    }


}
