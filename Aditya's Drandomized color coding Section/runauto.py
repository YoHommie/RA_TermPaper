import subprocess
import os

def compile_and_run_cpp(source_file, output_file="output.txt", executable_name="a.exe", iterations=20):
    # Compile the C++ file
    compile_command = ["g++", source_file, "-o", executable_name]
    
    try:
        subprocess.run(compile_command, check=True)
        print(f"Compilation successful: {executable_name}")
    except subprocess.CalledProcessError:
        print("Compilation failed!")
        return
    
    # Run the executable 20 times and capture output in append mode
    try:
        with open(output_file, "a") as out_file:
            for i in range(1, iterations + 1):
                out_file.write(f"\n--- Execution {i} ---\n")
                subprocess.run([executable_name], stdout=out_file, stderr=subprocess.STDOUT, check=True)
        print(f"All executions completed. Output saved to '{output_file}'.")
    except subprocess.CalledProcessError:
        print("Execution failed!")

if __name__ == "__main__":
    # Change this to your .cpp file as needed
    cpp_file = "P1T3_Using2stageHashing.cpp"
    output_file = "output.txt"
    executable_name = "a.exe"
    iterations = 20
    compile_and_run_cpp(cpp_file, output_file, executable_name, iterations)
