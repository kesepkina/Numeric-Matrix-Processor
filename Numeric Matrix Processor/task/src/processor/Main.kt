package processor

import java.lang.Exception
import java.math.RoundingMode
import kotlin.math.pow

enum class Operation(val menuNum: Int) {
    SUM(1),
    MULT_BY_CONST(2),
    MULT(3),
    TRANSPOSE(4),
    DETERMINANT(5),
    INVERSE(6),
    EXIT(0)
}

enum class TransposingType(val menuNum: Int) {
    MAIN_DIAG(1), SIDE_DIAG(2), VERTICAL_LINE(3), HORIZONTAL_LINE(4)
}

fun main() {
    while (true) {
        when (defineNextOperation()) {
            Operation.EXIT -> break
            Operation.SUM -> executeAddition()
            Operation.MULT -> executeMultiplication()
            Operation.MULT_BY_CONST -> executeMultiplicationByConstant()
            Operation.TRANSPOSE -> executeTransposing()
            Operation.DETERMINANT -> executeDeterminantCalculation()
            Operation.INVERSE -> executeInversion()
        }
        println()
    }


}

fun executeInversion() {
    val matrix = getMatrix()
    try {
        printMatrix(inverseMatrix(matrix))
    } catch (e: Exception) {
        println(e.message)
    }
}

fun inverseMatrix(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    if (!isSquare(matrix)) throw Exception("The operation cannot be performed.")
    val determinant = calcDeterminant(matrix)
    if (determinant == 0.0) throw Exception("This matrix doesn't have an inverse.")
    val adjoint = calcAdjoint(matrix)
    return mutiplyMatrixByNumber(executeTransposingByMainDiagonal(adjoint), 1 / determinant)
}

fun calcAdjoint(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    val resultMatrix = MutableList(matrix.size) { MutableList(matrix.size) { 0.0 } }
    for (i in resultMatrix.indices) {
        for (j in resultMatrix.indices) {
            resultMatrix[i][j] = calcMinor(matrix, i, j) * (-1.0).pow(i + j)
        }
    }
    return resultMatrix
}

fun executeDeterminantCalculation() {
    val matrix = getMatrix()
    if (!isSquare(matrix)) {
        println("The operation cannot be performed.")
        return
    }
    val determinant = calcDeterminant(matrix)
    println("The result is:\n$determinant")
}

fun calcDeterminant(matrix: MutableList<MutableList<Double>>): Double {
    if (matrix.size == 1) return matrix[0][0]
    if (matrix.size == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    var determinant = 0.0
    for (j in matrix[0].indices) {
        determinant += calcMinor(matrix, 0, j) * matrix[0][j] * (-1.0).pow(j)
    }
    return determinant
}

fun calcMinor(matrix: MutableList<MutableList<Double>>, i: Int, j: Int): Double {
    val minorMatrix = MutableList(matrix.size - 1) { MutableList(matrix.size - 1) { 0.0 } }
    var iMinor = 0
    for (ii in matrix.indices) {
        if (ii != i) {
            var jMinor = 0
            for (jj in matrix.indices) {
                if (jj != j) {
                    minorMatrix[iMinor][jMinor++] = matrix[ii][jj]
                }
            }
            iMinor++
        }
    }
    return calcDeterminant(minorMatrix)
}

fun executeTransposing() {
    val transpType = defineTransposingType()
    val matrix = getMatrix()
    printMatrix(
        when (transpType) {
            TransposingType.VERTICAL_LINE -> executeTransposingByVerticalLine(matrix)
            TransposingType.SIDE_DIAG -> executeTransposingBySideDiagonal(matrix)
            TransposingType.MAIN_DIAG -> executeTransposingByMainDiagonal(matrix)
            TransposingType.HORIZONTAL_LINE -> executeTransposingByHorizontalLine(matrix)
        }
    )
}

fun executeTransposingByHorizontalLine(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    val n = matrix.size
    val m = matrix[0].size
    val transposedMatrix = MutableList(n) { MutableList(m) { 0.0 } }
    for (i in transposedMatrix.indices) {
        for (j in transposedMatrix[0].indices) {
            transposedMatrix[i][j] = matrix[n - i - 1][j]
        }
    }
    return transposedMatrix
}

fun executeTransposingByVerticalLine(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    val n = matrix.size
    val m = matrix[0].size
    val transposedMatrix = MutableList(n) { MutableList(m) { 0.0 } }
    for (i in transposedMatrix.indices) {
        for (j in transposedMatrix[0].indices) {
            transposedMatrix[i][j] = matrix[i][m - j - 1]
        }
    }
    return transposedMatrix
}

fun executeTransposingBySideDiagonal(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    val n = matrix.size
    val m = matrix[0].size
    val transposedMatrix = MutableList(m) { MutableList(n) { 0.0 } }
    for (i in transposedMatrix.indices) {
        for (j in transposedMatrix[0].indices) {
            transposedMatrix[i][j] = matrix[m - j - 1][n - i - 1]
        }
    }
    return transposedMatrix
}

fun executeTransposingByMainDiagonal(matrix: MutableList<MutableList<Double>>): MutableList<MutableList<Double>> {
    val transposedMatrix = MutableList(matrix[0].size) { MutableList(matrix.size) { 0.0 } }
    for (i in transposedMatrix.indices) {
        for (j in transposedMatrix[0].indices) {
            transposedMatrix[i][j] = matrix[j][i]
        }
    }
    return transposedMatrix
}

fun defineTransposingType(): TransposingType {
    print(
        "\n1. Main diagonal\n" +
                "2. Side diagonal\n" +
                "3. Vertical line\n" +
                "4. Horizontal line\n" +
                "Your choice: "
    )
    return when (readln().toInt()) {
        TransposingType.HORIZONTAL_LINE.menuNum -> TransposingType.HORIZONTAL_LINE
        TransposingType.MAIN_DIAG.menuNum -> TransposingType.MAIN_DIAG
        TransposingType.SIDE_DIAG.menuNum -> TransposingType.SIDE_DIAG
        TransposingType.VERTICAL_LINE.menuNum -> TransposingType.VERTICAL_LINE
        else -> throw Exception("Undefined choice.")
    }
}

fun executeMultiplicationByConstant() {
    val matrix = getMatrix()
    val constant = getConstant()
    printMatrix(mutiplyMatrixByNumber(matrix, constant))
}

fun executeMultiplication() {
    val (matrix1, matrix2) = getTwoMatrices()
    try {
        printMatrix(calcMultOfMatrices(matrix1, matrix2))
    } catch (e: Exception) {
        println(e.message)
    }
}

fun executeAddition() {
    val (matrix1, matrix2) = getTwoMatrices()
    try {
        printMatrix(calcSumOfMatrices(matrix1, matrix2))
    } catch (e: Exception) {
        println(e.message)
    }
}

fun getMatrix(): MutableList<MutableList<Double>> {
    print("Enter size of matrix: ")
    val (n, m) = readln().split(" ").map { it.toInt() }
    println("Enter matrix:")
    val matrix = MutableList(n) { MutableList(m) { 0.0 } }
    for (i in 0 until n) {
        matrix[i] = readln().split(" ").map { it.toDouble() }.toMutableList()
    }
    return matrix
}

fun getConstant(): Double {
    print("Enter constant: ")
    val constant = readln().toDouble()
    return constant
}

fun getTwoMatrices(): Pair<MutableList<MutableList<Double>>, MutableList<MutableList<Double>>> {
    print("Enter size of first matrix: ")
    val (n1, m1) = readln().split(" ").map { it.toInt() }
    println("Enter first matrix:")
    val matrix1 = MutableList(n1) { MutableList(m1) { 0.0 } }
    for (i in 0 until n1) {
        matrix1[i] = readln().split(" ").map { it.toDouble() }.toMutableList()
    }
    print("Enter size of second matrix: ")
    val (n2, m2) = readln().split(" ").map { it.toInt() }
    println("Enter second matrix:")
    val matrix2 = MutableList(n2) { MutableList(m2) { 0.0 } }
    for (i in 0 until n2) {
        matrix2[i] = readln().split(" ").map { it.toDouble() }.toMutableList()
    }
    return Pair(matrix1, matrix2)
}

fun defineNextOperation(): Operation {
    print(
        "1. Add matrices\n" +
                "2. Multiply matrix by a constant\n" +
                "3. Multiply matrices\n" +
                "4. Transpose matrix\n" +
                "5. Calculate a determinant\n" +
                "6. Inverse matrix\n" +
                "0. Exit\n" +
                "Your choice: "
    )
    return when (readln().toInt()) {
        Operation.SUM.menuNum -> Operation.SUM
        Operation.MULT_BY_CONST.menuNum -> Operation.MULT_BY_CONST
        Operation.MULT.menuNum -> Operation.MULT
        Operation.EXIT.menuNum -> Operation.EXIT
        Operation.TRANSPOSE.menuNum -> Operation.TRANSPOSE
        Operation.DETERMINANT.menuNum -> Operation.DETERMINANT
        Operation.INVERSE.menuNum -> Operation.INVERSE
        else -> throw Exception("Your choice is undefined")
    }
}

fun mutiplyMatrixByNumber(matrix: MutableList<MutableList<Double>>, c: Double): MutableList<MutableList<Double>> {
    val resultMatrix = matrix.toMutableList()
    for (i in resultMatrix.indices) {
        for (j in resultMatrix[0].indices) {
            resultMatrix[i][j] = matrix[i][j] * c
        }
    }
    return resultMatrix
}

fun calcMultOfMatrices(
    matrix1: MutableList<MutableList<Double>>,
    matrix2: MutableList<MutableList<Double>>
): MutableList<MutableList<Double>> {
    val n1 = matrix1.size
    val m1 = matrix1[0].size
    val n2 = matrix2.size
    val m2 = matrix2[0].size
    if (!validateDimensionsForMult(n1, m1, n2, m2)) {
        throw Exception("The operation cannot be performed.")
    }
    val resultMatrix = MutableList(n1) { MutableList(m2) { 0.0 } }
    for (i in resultMatrix.indices) {
        for (j in resultMatrix[0].indices) {
            var element = 0.0
            for (k in matrix1[i].indices) {
                element += matrix1[i][k] * matrix2[k][j]
            }
            resultMatrix[i][j] = element
        }
    }
    return resultMatrix
}

fun calcSumOfMatrices(
    matrix1: MutableList<MutableList<Double>>,
    matrix2: MutableList<MutableList<Double>>
): MutableList<MutableList<Double>> {
    val n = matrix1.size
    val m = matrix1[0].size
    if (!validateDimensionsForSum(n, m, matrix2.size, matrix2[0].size)) {
        throw Exception("The operation cannot be performed.")
    }
    val sumMatrix = MutableList(n) { MutableList(m) { 0.0 } }
    for (i in sumMatrix.indices) {
        for (j in sumMatrix[0].indices) {
            sumMatrix[i][j] = matrix1[i][j] + matrix2[i][j]
        }
    }
    return sumMatrix
}

fun validateDimensionsForSum(rowsNum1: Int, colNum1: Int, rowsNum2: Int, colNum2: Int) =
    !(rowsNum1 != rowsNum2 || colNum1 != colNum2)

fun validateDimensionsForMult(rowsNum1: Int? = null, colNum1: Int, rowsNum2: Int, colNum2: Int? = null) =
    colNum1 == rowsNum2

fun isSquare(matrix: MutableList<MutableList<Double>>) = matrix.size == matrix[0].size

fun printMatrix(matrix: MutableList<MutableList<Double>>, prefix: String = "The result is:\n") {
    print(prefix)
    matrix.forEach { row ->
        println(
            row.joinToString(
                " ",
                transform = {
                    "%6s".format(
                        it.toBigDecimal().setScale(2, RoundingMode.HALF_DOWN).toString().trimEnd { last -> last == '0' }
                            .trimEnd { last -> last == '.' })
                })
        )
    }
}
