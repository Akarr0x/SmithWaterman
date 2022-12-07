#' Smith - Waterman
#'
#' Returns one of the optimal local aligments between two sequences
#' @param match is the match score in the smith waterman algorithm
#' @param missmatch is the missmatch score in the smith waterman algorithm
#' @param gap is the gap score in the smith waterman algorithm
#' @param first is the first sequence in the smith waterman algorithm
#' @param second is the second sequence in the smith waterman algorithm

#' @return The optimal alignment between the first sequence and the second,
#'         a result variable that contains a rappresentation of the differences
#' @examples 
#' optimal_alignment <- smithwaterman(DNAString("ACCTG"), DNAString("ATCTG"));
#' optimal_alignment <- smithwaterman(DNAString("ATCTG"), DNAString("ACCTG"), 3, -1, 0);
#' @export

smithwaterman <- function(first, second,
                          match = 3, missmatch = -3,
                          gap = -1) {
# Those are the first two sequences,
# while sub_matrix consists of the current scores:
#                               - Match score
#                               - Missmatch score
#                               - Gap penality
first_sequence <- first
second_sequence <- second
sub_matrix <- matrix(c(match, missmatch, gap), 1, 3)

if (length(first) == 0 || length(second) == 0)
    stop("Sequences length is not correct")
if (class(first_sequence) != "DNAString" ||
    class(second_sequence) != "DNAString")
    stop("Sequences must be DNAStrings")
if (!is.numeric(sub_matrix))
    stop("All the scores and penalities must be numeric")
if (match <= 0)
    stop("Match score must be higher than 0")
if (missmatch > 0)
    stop("Missmatch score must be lower than 0")
if (gap > 0)
    stop("Gap score must be lower than 0")

first_sequence <- unlist(strsplit(toString(first_sequence), ""))
second_sequence <- unlist(strsplit(toString(second_sequence), ""))
len1 <- length(first_sequence)
len2 <- length(second_sequence)
# This commands create the scoring matrix,
# that will be used to  compute the scores and calculate the maximum score
scoring_matrix <- matrix(, ncol = len2 + 1, nrow = len1 + 1)
rownames(scoring_matrix) <-  c("", unlist(strsplit(first_sequence, "")))
colnames(scoring_matrix) <-  c("", unlist(strsplit(second_sequence, "")))
# Initializing the scoring matrix with 0 as first row and first column,
# as the algorithm is designed, also initializing a traceback matrix that
# will be used later to correctly backtrace the maximum score
traceback_matrix <- scoring_matrix
scoring_matrix[, 1] <- 0
scoring_matrix[1, ] <- 0

for (i in 2:nrow(scoring_matrix)) {
    for (j in 2:ncol(scoring_matrix)){

        # The variable scores has the scores obtained from
        # an eventual match, a vertical and a horizzontal gap
        scores <-   c(scoring_matrix[i, j - 1] + sub_matrix[3],
                        scoring_matrix[i - 1, j] + sub_matrix[3])

        # The actual match is checked by this function
        if (colnames(scoring_matrix)[j] == rownames(scoring_matrix)[i])
            scores <- append(scores,
                        scoring_matrix[i - 1, j - 1] + sub_matrix[1])

        # In case the diagonal move does not have a match
        # when comparing the row and column names a missmatch must be added
        else
            scores <- append(scores,
                        scoring_matrix[i - 1, j - 1] + sub_matrix[2])

        # This function obtains the maximum from all the three scores
        scoring_matrix[i, j] <- max(scores)

        # As the algorithm suggests, if a value is < 0 we make it = 0
        if (scoring_matrix[i, j] <= 0)
                scoring_matrix[i, j] <- 0
    }
}

# I will now add a traceback function to
# correctly identify the different nucleotides matches
# Since multiple maximum scores are possible but all the different
# strings obtained (that end with this maximum) are considered equal
# in the case of multiple maximum score, I'll get the first one.
max_score <- which(scoring_matrix == max(scoring_matrix), arr.ind = TRUE)

if (length(max_score) > 2)
    max_score <- max_score[1, ]

traceback_score <- max(scoring_matrix)
traceback_position <- max_score
# Considering that the last value should always be a diagonal movement,
# since it's the only way to obtain a higher value, since othwerwhise we
# would have a gap penality
traceback_matrix[traceback_position[1], traceback_position[2]] <- 1

# This is the traceback function.
# It exits when one of two conditions are matched:
#   The traceback score is = 0, as said from the algorithm, if a 
#       score is = 0 the comparison is finished
#   The traceback hits position 1, in this case the first row and first
#       column consist in the "filler" column, that has only 0s
while (traceback_score != 0 ||
        traceback_position[1] != 1 && traceback_position[2] != 1) {
        largest <-
            c(scoring_matrix[traceback_position[1] - 1, traceback_position[2] - 1],
            scoring_matrix[traceback_position[1] - 1, traceback_position[2]],
            scoring_matrix[traceback_position[1], traceback_position[2] - 1])
        largest_ind <- which(largest == max(largest), arr.ind = TRUE)

        # In the case we have multiple maxes that means that we could
        # have possibly multiple comparisons with the same score.
        # Since it's needed only one possible alignment,
        # I will always take the first one. The way that the largest variable
        # is built, ensures that in the case of multiple paths, if
        # a diagonal movement is possible, it always chooses it.
        if (length(largest_ind) >= 2)
            largest_ind <- largest_ind[1]

        traceback_matrix[traceback_position[1],
                        traceback_position[2]] <- largest_ind

        # To develop a traceback I assumed that 1 means a diagonal movement
        # while 2 a row gap, and 3 a column gap
        if (largest_ind == 1) {
            traceback_score <- scoring_matrix[traceback_position[1] - 1,
                                             traceback_position[2] - 1]
            traceback_position[1] <- traceback_position[1] - 1
            traceback_position[2] <- traceback_position[2] - 1

        }
        if (largest_ind == 2) {
            traceback_score <- scoring_matrix[traceback_position[1] - 1,
                                             traceback_position[2]]
            traceback_position[1] <- traceback_position[1] - 1
        }
        if (largest_ind == 3) {
            traceback_score <- scoring_matrix[traceback_position[1],
                                             traceback_position[2] - 1]
            traceback_position[2] <- traceback_position[2] - 1
        }

        if (max(largest) == 0) {
            break
        }
}


# To actually obtain the optimal alignment between the two sequences
# a traceback algorithm for the traceback matrix is needed.
backtracking_score <- max_score

sequence_a <- c()
sequence_b <- c()

# This is the function used to backtrack the matrix.
# Starting from the start position, which is equal to the maximum score position
# it goes diagonally towards the origin point, stopping the first time it
# encounters a row, or a column, of which sum of values is equal to 0.
while (sum(!is.na(traceback_matrix[backtracking_score[1], ])) > 0 ||
        sum(!is.na(traceback_matrix[, backtracking_score[2]])) > 0) {

    # In case we have a sum of non NA values equals to 0 we do not have any
    # column gap, in this case we can add the column name to the first
    # sequence without add any gaps (that I described with '_')
    if (sum(!is.na(traceback_matrix[, backtracking_score[2]])) == 1) {
        sequence_a <- append(sequence_a, colnames(
                traceback_matrix)[backtracking_score[2]])
    }else {
        # In the case of a sum >= 2, that means that gaps are present,
        # so it's possible to add an amount of _ equal to the sum - 1,
        # since the last thing to add is the actual match, since the presence
        # of a gap indicates the presence of a match before it
        if (sum(!is.na(traceback_matrix[, backtracking_score[2]])) != 0) {
            sequence_a <- append(sequence_a, rep('_',
                (sum(!is.na(traceback_matrix[, backtracking_score[2]]))) - 1))

            sequence_a <- append(sequence_a, colnames(
                traceback_matrix)[backtracking_score[2]])
        }
    }
    if (sum(!is.na(traceback_matrix[backtracking_score[1], ])) == 1) {
        sequence_b <- append(sequence_b, rownames(
                traceback_matrix)[backtracking_score[1]])
    }else {
        if (sum(!is.na(traceback_matrix[backtracking_score[1], ])) != 0) {
            sequence_b <- append(sequence_b, rep('_', 
                (sum(!is.na(traceback_matrix[backtracking_score[1], ]))) - 1))

            sequence_b <- append(sequence_b, rownames(
                traceback_matrix)[backtracking_score[1]])
        }
    }
    if (backtracking_score[2] != 0)
        backtracking_score[2] <- backtracking_score[2] - 1
    if (backtracking_score[1] != 0)
        backtracking_score[1] <- backtracking_score[1] - 1

}

# It is also needed to reverse the results, since it
# starts comparing from the most far value

sequence_a <- rev(sequence_a)
if (sequence_a[1] == '')
    sequence_a <- sequence_a[-1]


sequence_b <- rev(sequence_b)
if (sequence_b[1] == '')
    sequence_b <- sequence_b[-1]


# This functions is needed to better show the differences between the two
# nucleotide strings, but could be omitted
# The visual rappresentation describes:
#            * In the case of a match between the two sequences
#            | In the case of a missmatch between the two sequences
#              In the case of a gap between the two sequences



results <- c()

for (value in seq_along(sequence_a)) {
    if (sequence_a[value] == sequence_b[value])
        results <- append(results, '*')
    else {
            if (sequence_a[value] == '_' || sequence_b[value] == '_')
                results <- append(results, ' ')
            else
                results <- append(results, '|')
        }
}
result_return <- list("first sequence" = sequence_b,
                      "second sequence" = sequence_a,
                      "result" = results, "score" = scoring_matrix[max_score])
return(result_return)
}
