import tkinter as tk
from tkinter.filedialog import askopenfilename
from collections import Counter
from pathlib import Path

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet Frequency")
        self.geometry("640x460")
        self.file_path = None

        
        btn_frame = tk.Frame(self)
        btn_frame.pack(pady=8)

        tk.Button(btn_frame, text="Pick File", width=20, command=self.pick_file).grid(row=0, column=0, padx=6)
        tk.Button(btn_frame, text="Process File", width=20, command=self.process_file).grid(row=0, column=1, padx=6)
        tk.Button(btn_frame, text="Reset", width=20, command=self.reset).grid(row=0, column=2, padx=6)

       
        self.status_label = tk.Label(self, text="No file selected.")
        self.status_label.pack(pady=6)

        self.output_box = tk.Text(self, height=18, width=80)
        self.output_box.pack(padx=10, pady=10)

    def pick_file(self):
        filename = askopenfilename(
            filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
        )
        if filename:
            self.file_path = Path(filename)
            self.status_label.config(text=f"Selected: {self.file_path}")
        else:
            self.file_path = None
            self.status_label.config(text="No file selected.")

    def reset(self):
        self.output_box.delete("1.0", tk.END)
        self.status_label.config(text="Cleared. Pick a file to begin.")

    def process_file(self):
        if not self.file_path:
            self.status_label.config(text="No file selected.")
            return

        try:
            seq = self._read_fasta_sequence(self.file_path)
            if not seq:
                self.status_label.config(text="File has no sequence lines.")
                return

          
            alph = list(dict.fromkeys(seq))
            counts = Counter(seq)
            n = len(seq)
            freqs = [counts[ch] / n for ch in alph]

           
            self.output_box.delete("1.0", tk.END)
            self.output_box.insert(tk.END, "Unique characters (alphabet) in order of appearance:\n")
            self.output_box.insert(tk.END, ", ".join(alph) + "\n\n")

            self.output_box.insert(tk.END, "Character  Count  Relative_Frequency\n")
            self.output_box.insert(tk.END, "-------------------------------------\n")
            for ch, f in zip(alph, freqs):
                self.output_box.insert(tk.END, f"{ch:>9}  {counts[ch]:>5}  {f:.6f}\n")

            self.status_label.config(text=f"Processed {n} symbols from: {self.file_path.name}")

        except Exception as e:
            self.status_label.config(text=f"Error: {e}")

    @staticmethod
    def _read_fasta_sequence(path: Path) -> str:
        
        seq_parts = []
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                seq_parts.append(line)
        return "".join(seq_parts).upper()

def main() -> int:
    App().mainloop()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
