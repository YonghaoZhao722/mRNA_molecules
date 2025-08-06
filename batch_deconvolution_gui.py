import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from skimage.io import imread, imsave
import numpy as np
import RedLionfishDeconv as rl
import threading
import time
from pathlib import Path


class DeconvolutionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Batch Deconvolution Tool")
        self.root.geometry("600x600")
        
        # Variables
        self.input_dir = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.psf_path = tk.StringVar()
        self.iterations = tk.IntVar(value=30)
        self.compute_method = tk.StringVar(value="gpu")
        self.processing = False
        self.estimated_time = tk.StringVar(value="Estimated time: --:--")
        
        self.create_widgets()
        
    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Input directory selection
        ttk.Label(main_frame, text="Input Directory:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.input_dir, width=50).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="Browse", command=self.select_input_dir).grid(row=0, column=2, padx=5, pady=5)
        
        # Output directory selection
        ttk.Label(main_frame, text="Output Directory:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.output_dir, width=50).grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="Browse", command=self.select_output_dir).grid(row=1, column=2, padx=5, pady=5)
        
        # PSF file selection
        ttk.Label(main_frame, text="PSF File:").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.psf_path, width=50).grid(row=2, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="Browse", command=self.select_psf_file).grid(row=2, column=2, padx=5, pady=5)
        
        # Iterations setting
        ttk.Label(main_frame, text="Iterations:").grid(row=3, column=0, sticky=tk.W, pady=5)
        iterations_frame = ttk.Frame(main_frame)
        iterations_frame.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Entry(iterations_frame, textvariable=self.iterations, width=10).pack(side=tk.LEFT)
        
        # Compute method setting (GPU/CPU)
        ttk.Label(main_frame, text="Compute Method:").grid(row=4, column=0, sticky=tk.W, pady=5)
        method_combo = ttk.Combobox(main_frame, textvariable=self.compute_method, width=15, state="readonly")
        method_combo['values'] = ('gpu', 'cpu')
        method_combo.grid(row=4, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Process button
        self.process_btn = ttk.Button(main_frame, text="Start Processing", command=self.start_processing)
        self.process_btn.grid(row=5, column=1, pady=20)
        
        # Estimated time label
        self.time_label = ttk.Label(main_frame, textvariable=self.estimated_time)
        self.time_label.grid(row=6, column=0, columnspan=3, pady=(5, 0))
        
        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='determinate')
        self.progress.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        # Status text
        self.status_text = tk.Text(main_frame, height=10, width=70)
        self.status_text.grid(row=8, column=0, columnspan=3, pady=10, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Scrollbar for status text
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=self.status_text.yview)
        scrollbar.grid(row=8, column=3, sticky=(tk.N, tk.S), pady=10)
        self.status_text.configure(yscrollcommand=scrollbar.set)
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(8, weight=1)
        
    def select_input_dir(self):
        directory = filedialog.askdirectory(title="Select Input Directory")
        if directory:
            self.input_dir.set(directory)
            self.log_message(f"Input directory selected: {directory}")
            
    def select_output_dir(self):
        directory = filedialog.askdirectory(title="Select Output Directory")
        if directory:
            self.output_dir.set(directory)
            self.log_message(f"Output directory selected: {directory}")
            
    def select_psf_file(self):
        file_path = filedialog.askopenfilename(
            title="Select PSF File",
            filetypes=[("TIFF files", "*.tif *.tiff"), ("All files", "*.*")]
        )
        if file_path:
            self.psf_path.set(file_path)
            self.log_message(f"PSF file selected: {file_path}")
            
    def log_message(self, message):
        self.status_text.insert(tk.END, message + "\n")
        
        # Limit text widget to prevent memory issues and improve performance
        max_lines = 100
        lines = self.status_text.get("1.0", tk.END).count('\n')
        if lines > max_lines:
            # Remove old lines from the top
            lines_to_remove = lines - max_lines
            self.status_text.delete("1.0", f"{lines_to_remove + 1}.0")
        
        self.status_text.see(tk.END)
        
        # Always update GUI for immediate feedback during processing
        self.root.update_idletasks()
        
    def validate_inputs(self):
        if not self.input_dir.get():
            messagebox.showerror("Error", "Please select an input directory")
            return False
            
        if not self.output_dir.get():
            messagebox.showerror("Error", "Please select an output directory")
            return False
            
        if not self.psf_path.get():
            messagebox.showerror("Error", "Please select a PSF file")
            return False
            
        if not os.path.exists(self.input_dir.get()):
            messagebox.showerror("Error", "Input directory does not exist")
            return False
            
        if not os.path.exists(self.psf_path.get()):
            messagebox.showerror("Error", "PSF file does not exist")
            return False
            
        if self.iterations.get() <= 0:
            messagebox.showerror("Error", "Iterations must be a positive number")
            return False
            
        return True
        
    def start_processing(self):
        if self.processing:
            return
            
        if not self.validate_inputs():
            return
            
        self.processing = True
        self.process_btn.config(text="Processing...", state="disabled")
        
        # Start processing in a separate thread to avoid freezing the GUI
        thread = threading.Thread(target=self.process_files)
        thread.daemon = True
        thread.start()
        
    def process_files(self):
        try:
            # Create output directory
            os.makedirs(self.output_dir.get(), exist_ok=True)
            self.log_message(f"Created output directory: {self.output_dir.get()}")
            
            # Load PSF
            self.log_message("Loading PSF...")
            psf = imread(self.psf_path.get())
            self.log_message("PSF loaded successfully")
            
            # Get list of TIFF files
            input_path = Path(self.input_dir.get())
            tiff_files = [f for f in input_path.iterdir() 
                         if f.suffix.lower() in ['.tif', '.tiff'] and not f.name.startswith('deconv_')]
            
            if not tiff_files:
                self.log_message("No TIFF files found in input directory")
                return
                
            self.log_message(f"Found {len(tiff_files)} TIFF files to process")
            
            # Setup progress bar
            self.progress.config(maximum=len(tiff_files))
            
            # Initialize timing variables
            start_time = time.time()
            
            # Process each file
            for i, file_path in enumerate(tiff_files):
                try:
                    # Log start of processing for each file
                    self.log_message(f"Processing: {file_path.name} ({i+1}/{len(tiff_files)})")
                    
                    # Load image
                    image = imread(str(file_path))
                    
                    # Perform deconvolution
                    compute_method = self.compute_method.get()
                    
                    deconvolved = rl.doRLDeconvolutionFromNpArrays(
                        image, psf, niter=self.iterations.get(), 
                        method=compute_method, resAsUint8=False
                    )
                    
                    # Keep result as float32 (default output)
                    result = deconvolved.astype(np.float32)
                    
                    # Save result
                    output_path = Path(self.output_dir.get()) / f'deconv_{file_path.name}'
                    imsave(str(output_path), result)
                    
                    # Log completion for each file
                    self.log_message(f"Completed: {output_path.name}")
                    
                except Exception as e:
                    self.log_message(f"Error processing {file_path.name}: {str(e)}")
                
                # Update progress bar after each file
                self.progress.config(value=i + 1)
                
                # Calculate and update estimated time after each file
                elapsed_time = time.time() - start_time
                avg_time_per_file = elapsed_time / (i + 1)
                remaining_files = len(tiff_files) - (i + 1)
                estimated_remaining = avg_time_per_file * remaining_files
                
                if remaining_files > 0:
                    hours = int(estimated_remaining // 3600)
                    minutes = int((estimated_remaining % 3600) // 60)
                    seconds = int(estimated_remaining % 60)
                    
                    if hours > 0:
                        time_str = f"Estimated time remaining: {hours:02d}:{minutes:02d}:{seconds:02d}"
                    else:
                        time_str = f"Estimated time remaining: {minutes:02d}:{seconds:02d}"
                else:
                    time_str = "Processing complete!"
                
                self.estimated_time.set(time_str)
                
            self.log_message("All files processed successfully!")
            messagebox.showinfo("Complete", "Batch deconvolution completed!")
            
        except Exception as e:
            self.log_message(f"Fatal error: {str(e)}")
            messagebox.showerror("Error", f"Processing failed: {str(e)}")
            
        finally:
            self.processing = False
            self.process_btn.config(text="Start Processing", state="normal")
            self.progress.config(value=0)
            self.estimated_time.set("Estimated time: --:--")


def main():
    root = tk.Tk()
    app = DeconvolutionGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()