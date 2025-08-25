import re


class Node:
    def __init__(self, value, parent=None):
        self.value = value
        self.children = []
        self.parent = parent

    def set_value(self, value):
        self.value = value

    # define a string representation of the node
    def __repr__(self):
        return self.value

class Newick:
    def __init__(self, input):
        self.root = Node(value="")
        self.breadcrumbs = [self.root]
        self.parse(input)

    def __repr__(self):
        result = "root: " + self.root.value
        result += "\nroot children: " + str(self.root.children)
        if len(self.root.children) > 0:
            result += "\nroot children children: " + str(self.root.children[0].children)
        return result

    def parse(self, input):
        before, character, after = self.scan(input)
        if character == "(":
            new_node = Node(value="", parent=self.breadcrumbs[-1])
            self.breadcrumbs[-1].children.append(new_node)
            self.breadcrumbs.append(new_node)
            self.parse(after)
        elif character == ",":
            if before:
                new_node = Node(value=before, parent=self.breadcrumbs[-1])
                self.breadcrumbs[-1].children.append(new_node)
            self.parse(after)
        elif character == ")":
            if before:
                new_node = Node(value=before, parent=self.breadcrumbs[-1])
                self.breadcrumbs[-1].children.append(new_node)

            # handle the case where a name appears after the )
            match = re.match(r"([A-Za-z0-9_]+)", after)
            if match:
                self.breadcrumbs[-1].value = match.group(1)
                after = after[match.end():]
            self.breadcrumbs.pop()

            self.parse(after)
        elif character == ";":
            if before:
                self.root.children[0].set_value(before)


    def scan(self, input):
        # Scan the input for the next special character and return the before, character, and after
        match = re.search(r'([(),;])', input)
        if match:
            before = input[:match.start()]
            character = match.group(1)
            after = input[match.end():]
            return before, character, after
        return input, None, None

    def search(self, value):
        # Search for a node with the given value
        return self._search(self.root, value)

    def _search(self, node, value):
        if node.value == value:
            return node
        for child in node.children:
            result = self._search(child, value)
            if result:
                return result
        return None
    
    def traverseChildren(self, start, target):
        # traverse the tree from start to target, counting edges
        if start.value == target:
            return 0
        for child in start.children:
            result = self.traverseChildren(child, target)
            if result is not None:
                return result + 1
        return None
    
    def dist(self, v1, v2):
        #print("searching for " + v1 + " to " + v2)
        v1_node = self.search(v1)
        #print("found v1 as " + str(v1_node))
        firstPass = self.traverseChildren(v1_node, v2)
        if firstPass is not None:
            #print(f"found on first pass! {firstPass}")
            return firstPass
        
        distance = 0
        current = v1_node
        matchFound = False
        while matchFound == False:
            distance += 1
            if current.parent is None:
                break
            current = current.parent
            if self.traverseChildren(current, v2) is not None:
                distance += self.traverseChildren(current, v2)
                matchFound = True
        return distance if matchFound else "search failed"


def main():
    
    dists = []

    # prompt for input
    file = input("Filename for input: ")
    file = re.sub(r'["\']', '', file)
    with open(file, "r") as f:
        lines = f.readlines()
        for i in range(0, len(lines), 3):
            data = lines[i:i+3]
            if data[0].endswith(";\n"):
                tree = Newick(data[0][:-1])
                #print(tree)
            if len(data) > 1 and len(data[1].strip().split()) == 2:
                v1, v2 = data[1].strip().split()
                dists.append(tree.dist(v1, v2))

    # print the distances with no newlines and spaces separating each value
    print(" ".join(map(str, dists)))
    
    '''
    input = "(dog,cat);"
    tree = Newick(input)
    print(tree)
    '''

if __name__ == "__main__":
    main()