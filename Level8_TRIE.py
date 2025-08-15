import re

class Tree:
    def __init__(self):
        self.numNodes = 1
        self.root = Node("root", 1)

    def add_sequence(self, s):
        current_node = self.root
        for char in s:
            if char not in current_node.children:
                current_node = current_node.add_child(char, self.numNodes + 1)
                self.numNodes += 1
            else:
                current_node = current_node.children[char]
            

    def print_tree(self, file = None):
        self.root.print_all_children(file)


class Node:
    def __init__(self, value, position):
        self.value = value
        self.children = {}
        self.position = position

    def add_child(self, char, position):
        if char not in self.children:
            self.children[char] = Node(char, position)
        return self.children[char]
    
    def print_all_children(self, file=None):
        for child in self.children.values():
            if file:
                file.write(f"{self.position} {child.position} {child.value}\n")
            else:
                print(f"{self.position} {child.position} {child.value}")
            child.print_all_children(file)

def main():

    # prompt for input
    file = input("Filename for input: ")
    file = re.sub(r'["\']', '', file)
    with open(file, "r") as f:
        data = f.read()

    seqs = data.strip().split()

    tree = Tree()
    for seq in seqs:
        tree.add_sequence(seq)

    # print tree into file in the same directory as the input file
    output_file = re.sub(r"[^\\]+\.txt$", 'output.txt', file)
    with open(output_file, "w") as f:
        tree.print_tree(f)


if __name__ == "__main__":
    main()