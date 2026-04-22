from typing import Any, Callable


def merge(l1: list, l2: list, *, key: Callable[[Any], Any] = lambda x: x) -> list:
    """Merge two lists preserving the relative order of elements from both inputs.

    Returns a list containing all elements from l1 and l2, deduplicated by key,
    ordered consistently with both inputs. Original items (not key values) are returned;
    when the same key appears in both lists, the item from l1 is kept.

    Args:
        l1: First ordered list.
        l2: Second ordered list.
        key: Function extracting a hashable identity from each element. Defaults to
             the element itself. Two elements with the same key are considered equal.

    Raises ValueError if the ordering constraints from l1 and l2 conflict,
    or if either input list contains duplicate elements (by key).
    """
    for lst, name in ((l1, "l1"), (l2, "l2")):
        keys = [key(item) for item in lst]
        if len(keys) != len(set(keys)):
            raise ValueError(f"Input {name} contains duplicate elements")

    # Assign a stable rank to each element keyed by key(item).
    # When the same key appears in both lists, the item from l1 is kept.
    rank: dict[Any, int] = {}
    items: dict[Any, Any] = {}  # key -> original item
    for item in l1 + l2:
        k = key(item)
        if k not in rank:
            rank[k] = len(rank)
            items[k] = item

    all_keys = list(rank.keys())
    graph: dict[Any, set] = {k: set() for k in all_keys}
    in_degree: dict[Any, int] = {k: 0 for k in all_keys}

    for lst in (l1, l2):
        for a, b in zip(lst, lst[1:]):
            ka, kb = key(a), key(b)
            if kb not in graph[ka]:
                graph[ka].add(kb)
                in_degree[kb] += 1

    queue = sorted([k for k in all_keys if in_degree[k] == 0], key=rank.__getitem__)
    result_keys: list = []

    while queue:
        k = queue.pop(0)
        result_keys.append(k)
        newly_free = []
        for successor in graph[k]:
            in_degree[successor] -= 1
            if in_degree[successor] == 0:
                newly_free.append(successor)
        queue = sorted(queue + newly_free, key=rank.__getitem__)

    if len(result_keys) < len(all_keys):
        remaining = [items[k] for k in all_keys if k not in set(result_keys)]
        raise ValueError(
            f"Cannot merge lists: ordering conflict detected among elements {remaining}"
        )
    return [items[k] for k in result_keys]


def merge_with(
    l1: list,
    l2: list,
    *,
    equiv: Callable[[Any, Any], bool],
    combine: Callable[[Any, Any], Any],
) -> list:
    """Merge two lists using a pairwise equivalence predicate and combine function.

    Like merge(), but instead of a key function, uses:
    - equiv(a, b): True if a (from l1) and b (from l2) should be treated as the same
      element and combined. Does not need to be an equivalence relation.
    - combine(a, b): produces the merged element when equiv is True.

    Matching is greedy: each element of l1 is matched to the first unmatched element
    of l2 that satisfies equiv. Unmatched elements pass through unchanged.
    No within-list duplicate checking is performed.

    Raises ValueError if ordering constraints from l1 and l2 conflict.
    """
    # Step 1: greedy matching — l1 drives, first-match wins
    matched_l1: dict[int, int] = {}  # l1 index -> l2 index
    matched_l2: dict[int, int] = {}  # l2 index -> l1 index
    available = list(range(len(l2)))
    for i, a in enumerate(l1):
        for pos, j in enumerate(available):
            if equiv(a, l2[j]):
                matched_l1[i] = j
                matched_l2[j] = i
                available.pop(pos)
                break

    # Step 2: build nodes in first-appearance order (l1 first, then unmatched l2)
    node_items: list[Any] = []
    l1_to_node: dict[int, int] = {}
    l2_to_node: dict[int, int] = {}

    for i, a in enumerate(l1):
        if i in matched_l1:
            j = matched_l1[i]
            idx = len(node_items)
            node_items.append(combine(a, l2[j]))
            l1_to_node[i] = idx
            l2_to_node[j] = idx
        else:
            idx = len(node_items)
            node_items.append(a)
            l1_to_node[i] = idx

    for j, b in enumerate(l2):
        if j not in matched_l2:
            idx = len(node_items)
            node_items.append(b)
            l2_to_node[j] = idx

    n = len(node_items)
    rank = {i: i for i in range(n)}
    graph: dict[int, set[int]] = {i: set() for i in range(n)}
    in_degree: dict[int, int] = {i: 0 for i in range(n)}

    # Step 3: ordering constraints from both lists
    for i in range(len(l1) - 1):
        a, b = l1_to_node[i], l1_to_node[i + 1]
        if a != b and b not in graph[a]:
            graph[a].add(b)
            in_degree[b] += 1

    for j in range(len(l2) - 1):
        a, b = l2_to_node[j], l2_to_node[j + 1]
        if a != b and b not in graph[a]:
            graph[a].add(b)
            in_degree[b] += 1

    # Step 4: topological sort, rank as tiebreaker
    queue = sorted([i for i in range(n) if in_degree[i] == 0], key=rank.__getitem__)
    result: list[int] = []

    while queue:
        node = queue.pop(0)
        result.append(node)
        newly_free = []
        for s in graph[node]:
            in_degree[s] -= 1
            if in_degree[s] == 0:
                newly_free.append(s)
        queue = sorted(queue + newly_free, key=rank.__getitem__)

    if len(result) < n:
        raise ValueError("Cannot merge lists: ordering conflict detected")

    return [node_items[i] for i in result]
