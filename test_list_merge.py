import pytest
from list_merge import merge, merge_with


@pytest.mark.parametrize("l1,l2,expected", [
    pytest.param(["b", "a"], ["a", "c"],           ["b", "a", "c"],     id="basic"),
    pytest.param(["x", "y"], ["z"],                ["x", "y", "z"],     id="unique_to_l1"),
    pytest.param(["z"],      ["x", "y"],           ["z", "x", "y"],     id="unique_to_l2"),
    pytest.param([],         [],                   [],                  id="both_empty"),
    pytest.param([],         ["a", "b"],           ["a", "b"],          id="l1_empty"),
    pytest.param(["a", "b"], [],                   ["a", "b"],          id="l2_empty"),
    pytest.param(["a", "b", "c"], ["a", "b", "c"], ["a", "b", "c"],    id="identical"),
    pytest.param(["a"],      ["b"],                ["a", "b"],          id="single_no_conflict"),
    pytest.param(["a"],      ["a"],                ["a"],               id="single_same"),
    pytest.param(["a", "c"], ["b", "c"],           ["a", "b", "c"],    id="ambiguous_deterministic"),
    pytest.param(["a", "b", "c"], ["a", "c"],      ["a", "b", "c"],    id="longer_chain"),
])
def test_merge(l1, l2, expected):
    assert merge(l1, l2) == expected


@pytest.mark.parametrize("l1,l2", [
    pytest.param(["a", "b"], ["b", "a"], id="conflict"),
])
def test_merge_raises_conflict(l1, l2):
    with pytest.raises(ValueError, match="conflict"):
        merge(l1, l2)


def test_within_list_duplicates():
    with pytest.raises(ValueError):
        merge(["a", "a"], ["a"])


# --- key= tests ---

def test_key_basic():
    # Objects compared by .id, not by equality
    class Item:
        def __init__(self, id, val):
            self.id = id
            self.val = val
        def __eq__(self, other):
            raise AssertionError("__eq__ must not be called")
        def __hash__(self):
            raise AssertionError("__hash__ must not be called")

    b = Item("b", 1)
    a1 = Item("a", 2)
    a2 = Item("a", 3)  # same key "a" as a1 — item from l1 (a1) should be kept
    c = Item("c", 4)

    result = merge([b, a1], [a2, c], key=lambda x: x.id)
    assert [x.id for x in result] == ["b", "a", "c"]
    assert result[1] is a1  # l1's item is kept when keys collide


def test_key_conflict():
    class Item:
        def __init__(self, id):
            self.id = id
        def __eq__(self, other):
            raise AssertionError("__eq__ must not be called")
        def __hash__(self):
            raise AssertionError("__hash__ must not be called")

    a, b = Item("a"), Item("b")
    with pytest.raises(ValueError, match="conflict"):
        merge([a, b], [b, a], key=lambda x: x.id)


def test_key_duplicate_in_input():
    class Item:
        def __init__(self, id):
            self.id = id

    a1, a2 = Item("a"), Item("a")
    with pytest.raises(ValueError):
        merge([a1, a2], [], key=lambda x: x.id)


# --- merge_with tests ---

def eq(a, b):
    return a == b

def keep(a, b):
    return a

def cat(a, b):
    return a + b


@pytest.mark.parametrize("l1,l2,expected", [
    # basic: no matches
    (["a", "b"], ["c", "d"],       ["a", "b", "c", "d"]),
    # single match in middle: combined, ordering preserved
    (["a", "b"], ["b", "c"],       ["a", "b", "c"]),
    # match at start
    (["a", "b"], ["a", "c"],       ["a", "b", "c"]),
    # match at end
    (["a", "b"], ["c", "b"],       ["a", "c", "b"]),
    # all match
    (["a", "b"], ["a", "b"],       ["a", "b"]),
    # empty lists
    ([],         ["a"],            ["a"]),
    (["a"],      [],               ["a"]),
    # combine is called: cat concatenates
    (["a"],      ["a"],            ["aa"]),
    # greedy: l1[0] matches l2[0], l1[1] gets no match (l2[0] already taken)
    (["a", "a"], ["a"],            ["a", "a"]),
])
def test_merge_with(l1, l2, expected):
    fn = cat if expected == ["aa"] else keep
    assert merge_with(l1, l2, equiv=eq, combine=fn) == expected


def test_merge_with_conflict():
    with pytest.raises(ValueError, match="conflict"):
        merge_with(["a", "b"], ["b", "a"], equiv=eq, combine=keep)


def test_merge_with_combine_called():
    combined = []
    def record_combine(a, b):
        combined.append((a, b))
        return f"{a}+{b}"
    result = merge_with(["x", "y"], ["y", "z"], equiv=eq, combine=record_combine)
    assert result == ["x", "y+y", "z"]
    assert combined == [("y", "y")]
