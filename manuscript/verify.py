"""Exact finite checks for the weighted Ev character-average identity."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial, prod
from math import gcd
from itertools import permutations


MAX_CHARACTER_N = 8
MAX_CYCLE_N = 60
MAX_Q_DEGREE = 10


@lru_cache(maxsize=None)
def partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None:
        largest = total
    return tuple(
        (first,) + tail
        for first in range(min(total, largest), 0, -1)
        for tail in partitions(total - first, first)
    )


def centralizer(partition):
    return prod(part ** count * factorial(count) for part, count in Counter(partition).items())


def ev_types(partition):
    choices = Counter({(): 1})
    for part in partition:
        updated = Counter()
        for cycle_type, multiplicity in choices.items():
            updated[tuple(sorted(cycle_type + (2 * part,), reverse=True))] += multiplicity
            updated[tuple(sorted(cycle_type + (part, part), reverse=True))] += multiplicity
        choices = updated
    assert sum(choices.values()) == 2 ** len(partition)
    return choices


@lru_cache(maxsize=None)
def rim_removals(shape, length):
    remaining_size = sum(shape) - length
    if remaining_size < 0:
        return ()
    removals = []
    for remainder in partitions(remaining_size):
        if len(remainder) > len(shape):
            continue
        padded = remainder + (0,) * (len(shape) - len(remainder))
        if any(smaller > larger for smaller, larger in zip(padded, shape)):
            continue
        cells = {(row, column) for row, width in enumerate(shape)
                 for column in range(padded[row], width)}
        if len(cells) != length:
            continue
        if any({(row, column), (row + 1, column), (row, column + 1),
                (row + 1, column + 1)} <= cells for row, column in cells):
            continue
        frontier = [next(iter(cells))]
        connected = set(frontier)
        while frontier:
            row, column = frontier.pop()
            for neighbor in [(row - 1, column), (row + 1, column),
                             (row, column - 1), (row, column + 1)]:
                if neighbor in cells and neighbor not in connected:
                    connected.add(neighbor)
                    frontier.append(neighbor)
        if connected == cells:
            height = len({row for row, column in cells}) - 1
            removals.append((remainder, (-1) ** height))
    return tuple(removals)


@lru_cache(maxsize=None)
def character(shape, cycle_type):
    assert sum(shape) == sum(cycle_type)
    if not cycle_type:
        return int(not shape)
    return sum(sign * character(remainder, cycle_type[1:])
               for remainder, sign in rim_removals(shape, cycle_type[0]))


def dimension(shape):
    hooks = [width - column + sum(other > column for other in shape[row + 1:])
             for row, width in enumerate(shape) for column in range(width)]
    return factorial(sum(shape)) // prod(hooks)


def laurent_product(partition):
    polynomial = {0: 1}
    for part in partition:
        updated = defaultdict(int)
        for exponent, coefficient in polynomial.items():
            for shift in [-part, 0, part]:
                updated[exponent + shift] += 2 * coefficient
        polynomial = dict(updated)
    return polynomial


def common_value(partition):
    polynomial = laurent_product(partition)
    return polynomial.get(0, 0) - polynomial.get(1, 0)


def closed_form(degree):
    return comb(degree // 2 + 2, 2) if degree % 2 == 0 else 0


def check_characters():
    root_pairs = set()
    partition_count = 0
    averages = []
    for degree in range(MAX_CHARACTER_N + 1):
        even_shapes = tuple(tuple(2 * part for part in shape)
                            for shape in partitions(degree) if len(shape) <= 3)
        rectangle = (degree, degree) if degree else ()
        for shape in set(even_shapes) | {rectangle}:
            assert character(shape, (1,) * (2 * degree)) == dimension(shape)
        left_average = Fraction(0)
        right_average = Fraction(0)
        for partition in partitions(degree):
            partition_count += 1
            left_sum = 0
            right_sum = 0
            for cycle_type, multiplicity in ev_types(partition).items():
                for shape in even_shapes:
                    root_pairs.add((shape, cycle_type))
                root_pairs.add((rectangle, cycle_type))
                left_sum += multiplicity * (-1) ** len(cycle_type) * sum(
                    character(shape, cycle_type) for shape in even_shapes)
                right_sum += multiplicity * character(rectangle, cycle_type)
            assert left_sum == right_sum == common_value(partition), (partition, left_sum, right_sum)
            left_average += Fraction(left_sum, centralizer(partition))
            right_average += Fraction(right_sum, centralizer(partition))
        assert left_average == right_average == closed_form(degree), degree
        averages.append(str(left_average))
    print(f"Murnaghan-Nakayama character sums: all {partition_count} partitions with 0 <= n <= {MAX_CHARACTER_N}.")
    print(f"Distinct shape/cycle-type pairs compared: {len(root_pairs)}.")
    print("Identity-class values agree with the hook-length formula for every selected shape.")
    print("Both weighted averages: " + ", ".join(averages) + ".")


def check_cycle_series():
    cycle = [{0: Fraction(1)}]
    for degree in range(1, MAX_CYCLE_N + 1):
        current = defaultdict(Fraction)
        for part in range(1, degree + 1):
            for exponent, coefficient in cycle[degree - part].items():
                for shift in [-part, 0, part]:
                    current[exponent + shift] += Fraction(2, degree) * coefficient
        cycle.append(dict(current))
    for degree, current in enumerate(cycle):
        geometric = defaultdict(int)
        for neutral in range(degree + 1):
            for positive in range(degree - neutral + 1):
                negative = degree - neutral - positive
                geometric[positive - negative] += (neutral + 1) * (positive + 1) * (negative + 1)
        assert current == geometric, degree
        assert current.get(0, 0) - current.get(1, 0) == closed_form(degree), degree
    print(f"Cycle-index recurrence versus three geometric factors: every Laurent coefficient for 0 <= n <= {MAX_CYCLE_N}.")
    print(f"Adjacent-coefficient cancellation and parity formula: 0 <= n <= {MAX_CYCLE_N}.")


def multiply_truncated(left, right, degree):
    product = [Fraction(0)] * (degree + 1)
    for first, left_coefficient in enumerate(left):
        for second, right_coefficient in enumerate(right[:degree + 1 - first]):
            product[first + second] += left_coefficient * right_coefficient
    return product


def check_q_series():
    coefficients = [Fraction(0)] * (MAX_Q_DEGREE + 1)
    for degree in range(MAX_Q_DEGREE + 1):
        for partition in partitions(degree):
            weight = [Fraction(1)] + [Fraction(0)] * MAX_Q_DEGREE
            for part in partition:
                factor = [Fraction(0)] * (MAX_Q_DEGREE + 1)
                for multiple in range(1, MAX_Q_DEGREE // part + 1):
                    factor[multiple * part] = (-1) ** (multiple - 1)
                weight = multiply_truncated(weight, factor, MAX_Q_DEGREE)
            value = Fraction(common_value(partition), centralizer(partition))
            coefficients = [old + value * term for old, term in zip(coefficients, weight)]
    expected = [1, 0, 3, -4, 9, -12, 22, -36, 60, -88, 135]
    assert coefficients == expected
    print("N=1 q-series, coefficients of q^0 through q^10: " + ", ".join(map(str, coefficients)) + ".")


def check_same_parity_detector():
    vandermonde = {
        (2, 1, 0): 1, (2, 0, 1): -1, (1, 2, 0): -1,
        (1, 0, 2): 1, (0, 2, 1): 1, (0, 1, 2): -1,
    }
    checked = 0
    for degree in range(19):
        for first in range(degree + 1):
            for second in range(min(first, degree - first) + 1):
                third = degree - first - second
                if not 0 <= third <= second:
                    continue
                monomials = set(permutations((first, second, third)))
                detector = 0
                for monomial in monomials:
                    for shift, sign in vandermonde.items():
                        shape = tuple(monomial[index] + shift[index] - (2 - index)
                                      for index in range(3))
                        if (shape[0] >= shape[1] >= shape[2] >= 0
                                and shape[0] % 2 == shape[1] % 2 == shape[2] % 2):
                            detector += sign
                principal = sum(int(2 * monomial[0] + monomial[1] == degree)
                                - int(2 * monomial[0] + monomial[1] == degree + 1)
                                for monomial in monomials)
                assert detector == principal, (degree, first, second, third)
                checked += 1
    print(f"Same-parity detector: all {checked} monomial-symmetric basis elements through degree 18, including odd degrees.")


def check_corollaries():
    paths = {0: 1}
    riordan = []
    for degree in range(31):
        riordan.append(paths.get(0, 0))
        tableau_sum = 0
        for first in range(degree + 1):
            for second in range(min(first, degree - first) + 1):
                third = degree - first - second
                if 0 <= third <= second and first % 2 == second % 2 == third % 2:
                    tableau_sum += dimension(tuple(part for part in (first, second, third) if part))
        assert tableau_sum == riordan[-1], degree
        updated = defaultdict(int)
        for height, count in paths.items():
            updated[height + 1] += count
            if height:
                updated[height] += count
                updated[height - 1] += count
        paths = dict(updated)
    for repetitions in range(1, 16):
        assert common_value((1,) * repetitions) == 2 ** repetitions * riordan[repetitions]
        trinomial = sum(comb(repetitions, 2 * chosen) * comb(2 * chosen, chosen)
                        for chosen in range(repetitions // 2 + 1))
        for part in range(2, 6):
            assert common_value((part,) * repetitions) == 2 ** repetitions * trinomial
    gcd_cases = 0
    for degree in range(1, 13):
        for partition in partitions(degree):
            if gcd(*partition) > 1:
                polynomial = laurent_product(partition)
                assert polynomial.get(1, 0) == 0
                assert common_value(partition) == polynomial[0] >= 0
                gcd_cases += 1
    row_cases = 0
    for degree in range(1, 6):
        for partition in partitions(degree):
            for shape in partitions(2 * degree):
                if shape[0] > degree:
                    total = sum(multiplicity * (-1) ** len(cycle_type)
                                * character(shape, cycle_type)
                                for cycle_type, multiplicity in ev_types(partition).items())
                    assert total == 0, (partition, shape, total)
                    row_cases += 1
    print("Riordan paths versus same-parity tableau degrees: 0 <= m <= 30.")
    print("Equal-part formulas: 1 <= repetitions <= 15; 1 <= part <= 5.")
    print(f"Common-divisor corollary: all {gcd_cases} applicable partitions of sizes 1 through 12.")
    print(f"First-row vanishing: all {row_cases} partition/shape pairs with 1 <= n <= 5 and first row > n.")
    cycle = [{0: Fraction(1)}]
    for degree in range(1, 11):
        current = defaultdict(Fraction)
        for part in range(1, degree + 1):
            weight = (-1) ** part * (part + 1)
            for exponent, coefficient in cycle[degree - part].items():
                for shift in (-part, 0, part):
                    current[exponent + shift] += Fraction(2 * weight, degree) * coefficient
        cycle.append(dict(current))
        direct = sum(Fraction(common_value(partition), centralizer(partition))
                     * prod((-1) ** part * (part + 1) for part in partition)
                     for partition in partitions(degree))
        assert direct == current.get(0, 0) - current.get(1, 0), degree
    print("Arbitrary-weight cycle identity: w_k=(-1)^k(k+1), through degree 10.")


def check_source_counterexample():
    partition = (5, 2, 1)
    degree = sum(partition)
    even_rows = [tuple(2 * part for part in shape)
                 for shape in partitions(degree) if len(shape) <= 7]
    paired_rows = [tuple(part for part in shape for repeat in range(2))
                   for shape in partitions(degree) if len(shape) <= 3]
    left = sum(multiplicity * (-1) ** len(cycle_type)
               * sum(character(shape, cycle_type) for shape in even_rows)
               for cycle_type, multiplicity in ev_types(partition).items())
    right = sum(multiplicity * sum(character(shape, cycle_type) for shape in paired_rows)
                for cycle_type, multiplicity in ev_types(partition).items())
    assert left - right == 8, (left, right)
    print(f"Source Appendix A, lambda=(5,2,1), N=3: left={left}, right={right}, difference=8.")


def check_rank_two():
    checked = 0
    for degree in range(7):
        even_rows = [tuple(2 * part for part in shape)
                     for shape in partitions(degree) if len(shape) <= 5]
        paired_rows = [tuple(part for part in shape for repeat in range(2))
                       for shape in partitions(degree) if len(shape) <= 2]
        for partition in partitions(degree):
            left = sum(multiplicity * (-1) ** len(cycle_type)
                       * sum(character(shape, cycle_type) for shape in even_rows)
                       for cycle_type, multiplicity in ev_types(partition).items())
            right = sum(multiplicity * sum(character(shape, cycle_type) for shape in paired_rows)
                        for cycle_type, multiplicity in ev_types(partition).items())
            assert left == right, (partition, left, right)
            checked += 1
    print(f"N=2 pointwise identity: all {checked} partitions with 0 <= n <= 6, by Murnaghan-Nakayama.")


def main():
    check_characters()
    check_cycle_series()
    check_q_series()
    check_same_parity_detector()
    check_corollaries()
    check_source_counterexample()
    check_rank_two()
    assert common_value((2, 1)) == -4
    print("Mixed-part example lambda=(2,1): both character sums equal -4.")
    print("All checks passed using integers and exact rational arithmetic.")
    print("The all-degree statements follow from the written proofs in note.pdf.")


if __name__ == '__main__':
    main()
