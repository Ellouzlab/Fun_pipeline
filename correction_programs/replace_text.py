import argparse

def arg_parser():
    parser = argparse.ArgumentParser(description="Replace text in a file")
    parser.add_argument("--input", help="File to modify")
    parser.add_argument("--output", help="File to write to")
    parser.add_argument("--find", help="Text to replace")
    parser.add_argument("--replacement", help="Replacement text")
    return parser.parse_args()

def read_file(file):
    with open(file, "r") as f:
        return f.read()

def replace_text(file, old, new):
    text=read_file(file)
    return text.replace(old, new)

def write_file(file, text):
    with open(file, "w") as f:
        f.write(text)

def main():
    args = arg_parser()
    text = replace_text(args.input, args.find, args.replacement)
    write_file(args.output, text)

if __name__ == '__main__':
    main()