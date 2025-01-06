class PriorityQueue {
    constructor() {
        this.values = [];
    }

    enqueue(val, priority) {
        this.values.push({ val, priority });
        this.sort();
    }

    dequeue() {
        return this.values.shift();
    }

    sort() {
        this.values.sort((a, b) => a.priority - b.priority);
    }
}

class WeightedGraph {
    constructor() {
        this.adjacencyList = {};
    }

    addVertex(vertex) {
        if (!this.adjacencyList[vertex]) this.adjacencyList[vertex] = [];
    }

    addEdge(vertex1, vertex2, weight) {
        this.adjacencyList[vertex1].push({ node: vertex2, weight });
        this.adjacencyList[vertex2].push({ node: vertex1, weight });
    }

    dijkstra(start, finish) {
        const nodes = new PriorityQueue();
        const distances = {};
        const previous = {};
        let path = []; // to return at end
        let smallest;

        // build up initial state
        for (let vertex in this.adjacencyList) {
            if (vertex === start) {
                distances[vertex] = 0;
                nodes.enqueue(vertex, 0);
            } else {
                distances[vertex] = Infinity;
                nodes.enqueue(vertex, Infinity);
            }
            previous[vertex] = null;
        }

        // as long as there is something to visit
        while (nodes.values.length) {
            smallest = nodes.dequeue().val;
            if (smallest === finish) {
                // we are done
                // build up path to return at end
                while (previous[smallest]) {
                    path.push(smallest);
                    smallest = previous[smallest];
                }
                break;
            }

            if (smallest || distances[smallest] !== Infinity) {
                for (let neighbor in this.adjacencyList[smallest]) {
                    // find neighboring node
                    let nextNode = this.adjacencyList[smallest][neighbor];
                    // calculate new distance to neighboring node
                    let candidate = distances[smallest] + nextNode.weight;
                    let nextNeighbor = nextNode.node;
                    if (candidate < distances[nextNeighbor]) {
                        // updating new smallest distance to neighbor
                        distances[nextNeighbor] = candidate;
                        // updating previous - How we got to neighbor
                        previous[nextNeighbor] = smallest;
                        // enqueue in priority queue with new priority
                        nodes.enqueue(nextNeighbor, candidate);
                    }
                }
            }
        }

        return path.concat(smallest).reverse();
    }
}

// Example usage:
const graph = new WeightedGraph();
graph.addVertex("1");
graph.addVertex("2");
graph.addVertex("3");
graph.addVertex("4");
graph.addVertex("5");
graph.addVertex("6");

graph.addEdge("1", "2", 4);
graph.addEdge("1", "3", 2);
graph.addEdge("2", "4", 7);
graph.addEdge("2", "5", 1);
graph.addEdge("3", "6", 5);
graph.addEdge("5", "6", 3);

// console.log(graph.dijkstra("1", "6"));


// ------------------------------------------------------------

function rotateMatrix(matrix) {
    const n = matrix.length;
    for (let layer = 0; layer < n / 2; layer++) {
        let first = layer;
        let last = n - 1 - layer;
        for (let i = first; i < last; i++) {
            let offset = i - first;
            // save top
            let top = matrix[first][i];
            // left -> top
            matrix[first][i] = matrix[last - offset][first];
            // bottom -> left
            matrix[last - offset][first] = matrix[last][last - offset];
            // right -> bottom
            matrix[last][last - offset] = matrix[i][last];
            // top -> right
            matrix[i][last] = top;
        }
    }
    return matrix;
}

// Example usage:
let matrix = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
];

console.log(rotateMatrix(matrix));


// -----------------------------------------------------------------------


class FreqStack {
    constructor() {
        this.freqMap = new Map();
        this.groupMap = new Map();
        this.maxFreq = 0;
    }

    push(val) {
        let freq = (this.freqMap.get(val) || 0) + 1;
        this.freqMap.set(val, freq);
        if (!this.groupMap.has(freq)) {
            this.groupMap.set(freq, []);
        }
        this.groupMap.get(freq).push(val);
        this.maxFreq = Math.max(this.maxFreq, freq);
    }

    pop() {
        let val = this.groupMap.get(this.maxFreq).pop();
        if (this.groupMap.get(this.maxFreq).length === 0) {
            this.maxFreq--;
        }
        this.freqMap.set(val, this.freqMap.get(val) - 1);
        return val;
    }
}

// Example usage:
const freqStack = new FreqStack();
freqStack.push(5);
freqStack.push(7);
freqStack.push(5);
freqStack.push(7);
freqStack.push(4);
freqStack.push(5);
console.log(freqStack.pop()); // 5
console.log(freqStack.pop()); // 7
console.log(freqStack.pop()); // 5
console.log(freqStack.pop()); // 4



// -----------------------------------------------------------------------

function maxSlidingWindow(nums, k) {
    let deque = [];
    let result = [];

    for (let i = 0; i < nums.length; i++) {
        while (deque.length && deque[0] < i - k + 1) {
            deque.shift();
        }
        while (deque.length && nums[deque[deque.length - 1]] < nums[i]) {
            deque.pop();
        }
        deque.push(i);
        if (i >= k - 1) {
            result.push(nums[deque[0]]);
        }
    }

    return result;
}

// Example usage:
let arr = [1, 3, -1, -3, 5, 3, 6, 7];
let k = 3;
console.log(maxSlidingWindow(arr, k)); // [3, 3, 5, 5, 6, 7]

// -----------------------------------------------------------------------

class ListNode {
    constructor(value) {
        this.value = value;
        this.next = null;
    }
}

function detectAndRemoveLoop(head) {
    let slow = head;
    let fast = head;

    while (fast !== null && fast.next !== null) {
        slow = slow.next;
        fast = fast.next.next;

        if (slow === fast) {
            removeLoop(slow, head);
            return;
        }
    }
}

function removeLoop(loopNode, head) {
    let ptr1 = head;
    let ptr2;

    // Find the start of the loop
    while (true) {
        ptr2 = loopNode;
        while (ptr2.next !== loopNode && ptr2.next !== ptr1) {
            ptr2 = ptr2.next;
        }

        if (ptr2.next === ptr1) {
            break;
        }

        ptr1 = ptr1.next;
    }

    // Remove the loop
    ptr2.next = null;
}

// Example usage:
let head = new ListNode(1);
head.next = new ListNode(2);
head.next.next = new ListNode(3);
head.next.next.next = new ListNode(4);
head.next.next.next.next = new ListNode(5);
head.next.next.next.next.next = head.next.next;

detectAndRemoveLoop(head);

// Print the modified linked list
let current = head;
while (current !== null) {
    console.log(current.value);
    current = current.next;
}

// -----------------------------------------------------------------------

class MedianFinder {
    constructor() {
        this.low = new PriorityQueue((a, b) => b - a); // Max-Heap
        this.high = new PriorityQueue((a, b) => a - b); // Min-Heap
    }

    addNum(num) {
        if (this.low.size() === 0 || num < this.low.peek()) {
            this.low.enqueue(num);
        } else {
            this.high.enqueue(num);
        }

        if (this.low.size() > this.high.size() + 1) {
            this.high.enqueue(this.low.dequeue());
        } else if (this.high.size() > this.low.size()) {
            this.low.enqueue(this.high.dequeue());
        }
    }

    findMedian() {
        if (this.low.size() > this.high.size()) {
            return this.low.peek();
        } else {
            return (this.low.peek() + this.high.peek()) / 2;
        }
    }
}

class PriorityQueue {
    constructor(compare) {
        this.values = [];
        this.compare = compare;
    }

    enqueue(val) {
        this.values.push(val);
        this.values.sort(this.compare);
    }

    dequeue() {
        return this.values.shift();
    }

    peek() {
        return this.values[0];
    }

    size() {
        return this.values.length;
    }
}

// Example usage:
const medianFinder = new MedianFinder();
medianFinder.addNum(6);
console.log(medianFinder.findMedian()); // 6
medianFinder.addNum(10);
console.log(medianFinder.findMedian()); // 8
medianFinder.addNum(2);
console.log(medianFinder.findMedian()); // 6
medianFinder.addNum(8);
console.log(medianFinder.findMedian()); // 7
medianFinder.addNum(4);
console.log(medianFinder.findMedian()); // 6
medianFinder.addNum(12);
console.log(medianFinder.findMedian()); // 7


// -----------------------------------------------------------------------


class TreeNode {
    constructor(value) {
        this.value = value;
        this.left = null;
        this.right = null;
    }
}

function lowestCommonAncestor(root, p, q) {
    if (!root || root === p || root === q) return root;
    let left = lowestCommonAncestor(root.left, p, q);
    let right = lowestCommonAncestor(root.right, p, q);
    if (left && right) return root;
    return left ? left : right;
}

// Example usage:
let root = new TreeNode(3);
root.left = new TreeNode(5);
root.right = new TreeNode(1);
root.left.left = new TreeNode(6);
root.left.right = new TreeNode(2);
root.right.left = new TreeNode(0);
root.right.right = new TreeNode(8);
root.left.right.left = new TreeNode(7);
root.left.right.right = new TreeNode(4);

let p = root.left; // Node 5
let q = root.right; // Node 1

console.log("LCA:", lowestCommonAncestor(root, p, q).value); // LCA: 3


// -----------------------------------------------------------------------

class Graph {
    constructor() {
        this.adjacencyList = {};
    }

    addVertex(vertex) {
        if (!this.adjacencyList[vertex]) this.adjacencyList[vertex] = [];
    }

    addEdge(vertex1, vertex2) {
        this.adjacencyList[vertex1].push(vertex2);
    }

    topologicalSortDFS() {
        const stack = [];
        const visited = new Set();

        const dfs = (vertex) => {
            if (visited.has(vertex)) return;
            visited.add(vertex);
            for (let neighbor of this.adjacencyList[vertex]) {
                dfs(neighbor);
            }
            stack.push(vertex);
        };

        for (let vertex in this.adjacencyList) {
            if (!visited.has(vertex)) {
                dfs(vertex);
            }
        }

        return stack.reverse();
    }

    topologicalSortKahn() {
        const inDegree = {};
        const queue = [];
        const result = [];

        for (let vertex in this.adjacencyList) {
            inDegree[vertex] = 0;
        }

        for (let vertex in this.adjacencyList) {
            for (let neighbor of this.adjacencyList[vertex]) {
                inDegree[neighbor]++;
            }
        }

        for (let vertex in inDegree) {
            if (inDegree[vertex] === 0) {
                queue.push(vertex);
            }
        }

        while (queue.length) {
            let vertex = queue.shift();
            result.push(vertex);

            for (let neighbor of this.adjacencyList[vertex]) {
                inDegree[neighbor]--;
                if (inDegree[neighbor] === 0) {
                    queue.push(neighbor);
                }
            }
        }

        return result;
    }
}

// Example usage:
const dag = new Graph();
dag.addVertex("1");
dag.addVertex("2");
dag.addVertex("3");
dag.addVertex("4");

dag.addEdge("1", "2");
dag.addEdge("1", "3");
dag.addEdge("3", "4");
dag.addEdge("2", "4");

console.log("Topological Order (DFS):", dag.topologicalSortDFS()); // [1, 3, 2, 4] or [1, 2, 3, 4]
console.log("Topological Order (Kahn's Algorithm):", dag.topologicalSortKahn()); // [1, 3, 2, 4] or [1, 2, 3, 4]


// -----------------------------------------------------------------------

function countDistinctSubsequences(s) {
    const MOD = 1e9 + 7;
    const n = s.length;
    const dp = new Array(n + 1).fill(0);
    dp[0] = 1;

    const lastOccurrence = {};

    for (let i = 1; i <= n; i++) {
        dp[i] = (dp[i - 1] * 2) % MOD;
        const char = s[i - 1];
        if (lastOccurrence[char] !== undefined) {
            dp[i] = (dp[i] - dp[lastOccurrence[char] - 1] + MOD) % MOD;
        }
        lastOccurrence[char] = i;
    }

    return dp[n];
}

// Example usage:
const string = "abcbac";
console.log("Number of Distinct Subsequences:", countDistinctSubsequences(string)); // 32