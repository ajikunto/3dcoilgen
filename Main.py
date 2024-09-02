# Created by KuntoAji
import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import random
import numpy as np
from matplotlib.patches import Circle
import ezdxf

class HelixApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("3D Coil Spring Generator V1.0 - R&D PT. APM Armada Suspension")

        # Mendapatkan ukuran layar
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        self.geometry(f"{int(screen_width)}x{int(screen_height)}+10+10")

        # Set theme and color
        ctk.set_appearance_mode("System")
        ctk.set_default_color_theme("blue")

        self.helix_data = []
        self.last_endpoint = {'x': 0, 'y': 0, 'z': 0, 'theta': 0}

        # Set up the input fields with default values and specific validation
        self.radius1_var = ctk.DoubleVar(value=30.0)
        self.radius2_var = ctk.DoubleVar(value=30.0)
        self.turns_var = ctk.DoubleVar(value=1.00)
        self.coil_height_var = ctk.IntVar(value=10)
        self.line_width_var = ctk.DoubleVar(value=10.0)

        # Default warna helix
        self.helix_color = "blue"

        # Konfigurasi grid utama untuk responsivitas
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        self.create_widgets()
        self.create_tabs()

        # Bind event untuk tombol ENTER
        self.bind("<Return>", lambda event: self.add_helix())

        # Tambahkan event handler untuk ketika jendela ditutup
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
    
    def create_widgets(self):
        # Input Frame
        input_frame = ctk.CTkFrame(self)
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # Configure column weights
        input_frame.columnconfigure(0, weight=1)
        input_frame.columnconfigure(1, weight=1)
        input_frame.columnconfigure(2, weight=1)
        input_frame.rowconfigure(6, weight=1)

        # Input fields with validation and labels
        ctk.CTkLabel(input_frame, text="Mean Radius Start (mm):").grid(row=0, column=0, padx=5, sticky='w')
        ctk.CTkEntry(input_frame, textvariable=self.radius1_var, width=80).grid(row=0, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Mean Radius End (mm):").grid(row=1, column=0, padx=5, sticky='w')
        ctk.CTkEntry(input_frame, textvariable=self.radius2_var, width=80).grid(row=1, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Turns (coil):").grid(row=2, column=0, padx=5, sticky='w')
        ctk.CTkEntry(input_frame, textvariable=self.turns_var, width=80).grid(row=2, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Pitch (mm):").grid(row=3, column=0, padx=5, sticky='w')
        ctk.CTkEntry(input_frame, textvariable=self.coil_height_var, width=80).grid(row=3, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Bar Diameter (mm):").grid(row=4, column=0, padx=5, sticky='w')
        ctk.CTkEntry(input_frame, textvariable=self.line_width_var, width=80).grid(row=4, column=1, padx=5, pady=5)

        # Radio button untuk warna sama
        self.color_choice = ctk.StringVar(value="same")
        ctk.CTkRadioButton(input_frame, text="Same Color", variable=self.color_choice, value="same").grid(row=5, column=0, padx=5, pady=5)

        # Radio button untuk warna random
        ctk.CTkRadioButton(input_frame, text="Random Color", variable=self.color_choice, value="random").grid(row=5, column=1, padx=5, pady=5)

        # Buttons
        ctk.CTkButton(input_frame, text="Add Coil", command=self.add_helix).grid(row=0, column=2, padx=6, pady=5, sticky="ew")
        ctk.CTkButton(input_frame, text="Reset", command=self.reset_helix).grid(row=1, column=2, padx=6, pady=5, sticky="ew")
        ctk.CTkButton(input_frame, text="Export 2D DXF", command=self.export_dxf).grid(row=2, column=2, padx=6, pady=5, sticky="ew")
        # Copy Button
        ctk.CTkButton(input_frame, text="Copy Table Data", command=self.copy_table_data).grid(row=3, column=2, padx=5, pady=5, sticky="ew")

        # Table
        self.table = ttk.Treeview(input_frame, columns=("Turns", "Coil Height", "Angel", "Acc Turns", "Acc Height", "Inner Dia", "Outer Dia"), show='headings', height=8)
        self.table.grid(row=6, column=0, columnspan=3, padx=5, pady=5, sticky="nsew")

        # Define columns
        self.table.heading("Turns", text="Turns")
        self.table.heading("Coil Height", text="Pitch")
        self.table.heading("Angel", text="Angel")
        self.table.heading("Acc Turns", text="Accu Turns")
        self.table.heading("Acc Height", text="Accu Height")
        self.table.heading("Inner Dia", text="Inner Dia")
        self.table.heading("Outer Dia", text="Outer Dia")

        # Set column widths and alignment
        self.table.column("Turns", anchor="center", width=40)
        self.table.column("Coil Height", anchor="center", width=50)
        self.table.column("Angel", anchor="center", width=40)
        self.table.column("Acc Turns", anchor="center", width=70)
        self.table.column("Acc Height", anchor="center", width=70)
        self.table.column("Inner Dia", anchor="center", width=70)
        self.table.column("Outer Dia", anchor="center", width=70)
        
    def create_tabs(self):
        # Frame untuk tab
        plot_frame = ctk.CTkFrame(self)
        plot_frame.grid(row=0, column=1, rowspan=2, padx=10, pady=10, sticky="nsew")

        # Konfigurasi plot_frame untuk memperluas sesuai dengan ukuran window
        plot_frame.grid_rowconfigure(0, weight=1)
        plot_frame.grid_columnconfigure(0, weight=1)

        # Buat notebook untuk tab panel
        self.notebook = ttk.Notebook(plot_frame)
        self.notebook.pack(expand=1, fill='both')

        # Tab untuk 2D plot
        self.tab_2d = ctk.CTkFrame(self.notebook)
        self.notebook.add(self.tab_2d, text="2D Plot")

        # Tab untuk 3D plot
        self.tab_3d = ctk.CTkFrame(self.notebook)
        self.notebook.add(self.tab_3d, text="3D Plot")

        # Buat figure untuk plot 3D
        self.fig_3d = Figure(figsize=(9, 9), constrained_layout=True)
        self.ax_3d = self.fig_3d.add_subplot(1, 1, 1, projection='3d')
        self.ax_3d.set_aspect('equal')
        
        self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, master=self.tab_3d)
        self.canvas_3d.get_tk_widget().pack(expand=1, fill='both')

        # Buat figure untuk plot 2D
        self.fig_2d = Figure(constrained_layout=True)
        self.ax_front = self.fig_2d.add_subplot(1, 2, 1)
        self.ax_front.set_aspect('equal')
        self.ax_2d = self.fig_2d.add_subplot(1, 2, 2)
        self.ax_2d.set_aspect('equal')

        self.canvas_2d = FigureCanvasTkAgg(self.fig_2d, master=self.tab_2d)
        self.canvas_2d.get_tk_widget().pack(expand=1, fill='both')
      
        # Inisialisasi untuk plot tampak depan
        self.ax_front.set_xlabel("X")
        self.ax_front.set_ylabel("Y")

        # Note di area figure   
        self.fig_2d.text(0.05, 0.95, "Note:", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')
        self.fig_2d.text(0.05, 0.90, "Hold Left click to ⟷ Pan", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')
        self.fig_2d.text(0.05, 0.85, "Mouse Scroll to ▲ Zoom in / ▼ Zoom out", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')

        # Note di area figure   
        self.fig_3d.text(0.05, 0.95, "Note:", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')
        self.fig_3d.text(0.05, 0.90, "Hold Left click to ↻ Rotate", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')
        self.fig_3d.text(0.05, 0.85, "Hold Right click to ▲ Zoom in / ▼ Zoom out", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')
        self.fig_3d.text(0.05, 0.80, "Hold Middle button to ⟷ Pan", transform=self.fig_2d.transFigure, fontsize=10, color='blue', ha='left')

        # Contoh plot 3D (bisa diganti dengan data yang sesuai)
        self.ax_3d.set_xlabel("X")
        self.ax_3d.set_ylabel("Y")
        self.ax_3d.set_zlabel("Z")
        self.canvas_3d.draw()

        # Tambahkan interaksi zoom in/out dan pan pada plot 2D
        self.canvas_2d.mpl_connect("scroll_event", self.on_scroll)
        self.canvas_2d.mpl_connect("button_press_event", self.on_press)
        self.canvas_2d.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas_2d.mpl_connect("button_release_event", self.on_release)
        self.canvas_2d.draw()

        # Variabel untuk panning
        self.pan_active = False
        self.last_press = None

    def on_scroll(self, event):
        """Zoom in/out berdasarkan event scroll"""
        base_scale = 1.2
        self.zoom(event, self.ax_front, base_scale)
        self.zoom(event, self.ax_2d, base_scale)

        # Ubah kursor saat zoom aktif
        if event.button == 'up':
            self.canvas_2d.get_tk_widget().config(cursor="sizing")  # Zoom in
        elif event.button == 'down':
            self.canvas_2d.get_tk_widget().config(cursor="sizing")  # Zoom out

        self.canvas_2d.draw()

    def zoom(self, event, ax, base_scale):
        curr_xlim = ax.get_xlim()
        curr_ylim = ax.get_ylim()

        xdata = event.xdata  # posisi mouse di data coords
        ydata = event.ydata

        if event.button == 'up':
            # Zoom in
            scale_factor = 1/base_scale
        elif event.button == 'down':
            # Zoom out
            scale_factor = base_scale
        else:
            return

        new_width = (curr_xlim[1] - curr_xlim[0]) * scale_factor
        new_height = (curr_ylim[1] - curr_ylim[0]) * scale_factor

        relx = (curr_xlim[1] - xdata) / (curr_xlim[1] - curr_xlim[0])
        rely = (curr_ylim[1] - ydata) / (curr_ylim[1] - curr_ylim[0])

        # Set new limits
        ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * relx])
        ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * rely])

    def on_press(self, event):
        """Aktifkan panning ketika tombol mouse ditekan"""
        if event.button == 1:
            self.pan_active = True
            self.last_press = event
            self.canvas_2d.get_tk_widget().config(cursor="fleur")

    def on_motion(self, event):
        """Panning berdasarkan gerakan mouse"""
        if not self.pan_active or event.inaxes is None or self.last_press is None:
            return

        dx = event.xdata - self.last_press.xdata
        dy = event.ydata - self.last_press.ydata

        ax = event.inaxes
        ax.set_xlim(ax.get_xlim() - dx)
        ax.set_ylim(ax.get_ylim() - dy)

        self.last_press = event
        self.canvas_2d.draw()

    def on_release(self, event):
        """Nonaktifkan panning ketika tombol mouse dilepas"""
        self.pan_active = False
        self.last_press = None
        self.canvas_2d.get_tk_widget().config(cursor="arrow")


    def add_helix(self):
        try:
            radius1 = self.radius1_var.get()
            radius2 = self.radius2_var.get()
            turns = self.turns_var.get()
            coil_height = self.coil_height_var.get()
            line_width = self.line_width_var.get()
        except (tk.TclError, ValueError):
            messagebox.showerror("Input Error", "Please enter valid numeric values.")
            return

        if turns <= 0 or radius1 <= 0 or radius2 <= 0 or line_width <= 0:
            messagebox.showerror("Input Error", "Turn or Radius or Bar diameter must be non zero.")
            return

        # Determine helix color
        if self.color_choice.get() == "random":
            self.helix_color = (random.random(), random.random(), random.random())
        else:
            self.helix_color = "blue"

        start_theta = self.last_endpoint['theta']
        start_z = self.last_endpoint['z']
        t = np.linspace(start_theta, start_theta + turns * 2 * np.pi, 1000)
        z = np.linspace(start_z, start_z + coil_height * turns, 1000)
        r = np.linspace(radius1, radius2, 1000)
        x = r * np.cos(t)
        y = r * np.sin(t)

        # Plot 3D
        self.ax_3d.plot(x, y, z, color=self.helix_color, linewidth=line_width)
        
        # Hilangkan grid dan warna plane
        self.ax_3d.xaxis.pane.fill = False  # Plane X menjadi transparan
        self.ax_3d.yaxis.pane.fill = False  # Plane Y menjadi transparan
        self.ax_3d.zaxis.pane.fill = False  # Plane Z menjadi transparan

        # Hitung jarak tiap axis
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        z_range = z.max() - z.min()

        # Cari jarak terbesar
        max_range = np.max([x_range, y_range, z_range])

        # Tentukan tengah dari tiap axis
        x_mid = (x.max() + x.min()) * 0.5
        y_mid = (y.max() + y.min()) * 0.5
        z_mid = (z.max() + z.min()) * 0.5

        # Set limits supaya aspect ratio tetap benar
        self.ax_3d.set_xlim(x_mid - max_range/2, x_mid + max_range/2)
        self.ax_3d.set_ylim(y_mid - max_range/2, y_mid + max_range/2)
        self.ax_3d.set_zlim(z_mid - max_range/2, z_mid + max_range/2)

        # Optional: Hilangkan grid lines
        self.ax_3d.grid(False)

        self.canvas_3d.draw()
        self.fig_3d.tight_layout()

        # Plot 2D
        self.ax_2d.plot(z, x, linewidth=1)  # Z-axis horizontal (rotated 90 degrees)
        for i in range(0, len(z), len(z)):  # Add circles at each turn
            circle = Circle((z[i], x[i]), radius=line_width/2, color='r', fill=False)
            self.ax_2d.add_patch(circle)

        self.canvas_2d.draw()
        self.fig_2d.tight_layout()

        # Proyeksi tampak depan
        self.ax_front.plot(y, x, linewidth=1)  # Plot tampak depan (Z view)
        for i in range(0, len(y), len(y)):  # Add circles at each turn
            circle = Circle((y[i], x[i]), radius=line_width/2, color='b', fill=False)
            self.ax_front.add_patch(circle)

        # Sinkronisasi rentang sumbu (xlim dan ylim)
        self.ax_front.set_xlim([x.min()-line_width, x.max() + line_width])
        self.ax_2d.set_xlim([x.min()-line_width, x.max() + line_width])
        
        y_min, y_max = min(y), max(y)
        self.ax_front.set_ylim([y_min-line_width, y_max + line_width])
        self.ax_2d.set_ylim([y_min-line_width, y_max + line_width])


        # Simpan data untuk helix
        self.helix_data.append((x, y, z, radius1, radius2, turns, coil_height, line_width, t))

        # Update nilai terakhir
        self.last_endpoint['x'] = x[-1]
        self.last_endpoint['y'] = y[-1]
        self.last_endpoint['z'] = z[-1]
        self.last_endpoint['theta'] = t[-1]

        # mencari suduk slope helix
        # Hitung keliling rata-rata helix
        mean_circumference = np.pi * (radius1 + radius2) / 2

        # Hitung sudut kemiringan helix dalam radian
        theta_slope = np.arctan(coil_height / mean_circumference)

        # Jika ingin dalam derajat:
        theta_slope_deg = round(np.degrees(theta_slope), 2)

        # Tambahkan data ke tabel
        self.table.insert("", "end", values=(
            turns, 
            coil_height, 
            theta_slope_deg,  # Menampilkan nilai theta tengah
            sum([d[5] for d in self.helix_data]), 
            sum([self.helix_data[i][5] * (self.helix_data[i][6] if i == 0 else self.helix_data[i][6]) for i in range(len(self.helix_data))]), 
            (radius1 + radius2) - (line_width),
            (radius1 + radius2) + (line_width)
        ))

        # Update radius1 dengan radius2 sebelumnya secara otomatis
        self.radius1_var.set(radius2)

    def reset_helix(self):
        self.ax_3d.cla()  # Clear the 3D plot
        self.ax_2d.cla()  # Clear the 2D plot
        self.ax_front.cla()  # Clear the 2D Front plot
        self.canvas_3d.draw()
        self.canvas_2d.draw()
        self.helix_data.clear()  # Clear the helix data
        self.table.delete(*self.table.get_children())  # Clear the table
        self.last_endpoint = {'x': 0, 'y': 0, 'z': 0, 'theta': 0}  # Reset last endpoint
        self.radius1_var.set(30.0)  # Reset radius1 to default value
        self.radius2_var.set(30.0)  # Reset radius2 to default value
        self.turns_var.set(1.00)
        self.coil_height_var.set(10)
        self.line_width_var.set(10.0)

    def export_dxf(self):
        # Ask the user to choose a location to save the file
        file_path = filedialog.asksaveasfilename(defaultextension=".dxf",
                                                filetypes=[("DXF files", "*.dxf")],
                                                title="Save DXF File")

        if not file_path:
            return  # If the user cancels the dialog, do nothing
        try:
            doc = ezdxf.new()
            msp = doc.modelspace()

            # Ekstrak data dari grafik 2D (Tampak depan) dan tambahkan ke DXF
            for line in self.ax_front.get_lines():
                y_data = line.get_xdata()
                x_data = line.get_ydata()
                # Tambahkan garis-garis helix ke file DXF
                for i in range(len(y_data) - 1):
                    msp.add_line((y_data[i], x_data[i]), (y_data[i + 1], x_data[i + 1]))

            # Ekstrak data dari grafik 2D (X-Z) dan tambahkan ke DXF
            for line in self.ax_2d.get_lines():
                x_data = line.get_xdata()
                z_data = line.get_ydata()
                # Rotasi 90 derajat: (x, z) -> (z, -x)
                rotated_x_data = x_data + (self.radius2_var.get() * 3)
                rotated_z_data = z_data
                # Tambahkan garis-garis helix ke file DXF
                for i in range(len(rotated_x_data) - 1):
                    msp.add_line((rotated_x_data[i], rotated_z_data[i]), (rotated_x_data[i + 1], rotated_z_data[i + 1]))

            # Tambahkan lingkaran pada setiap akhir turn
            for line in self.ax_2d.get_lines():
                x_data = line.get_xdata()
                z_data = line.get_ydata()
                # Rotasi 90 derajat: (x, z) -> (z, -x)
                rotated_x_data = x_data + (self.radius2_var.get() * 3)
                rotated_z_data = z_data

                indices = np.linspace(0, len(rotated_x_data) - 1, num=int(self.turns_var.get() * 2) + 1, dtype=int)
                for idx in indices:
                    msp.add_circle((rotated_x_data[idx], rotated_z_data[idx]), self.line_width_var.get() / 2)

            doc.saveas(file_path)
            messagebox.showinfo("Export Successful", f"DXF file has been exported successfully to {file_path}!")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Failed to export DXF: {e}")
    
    def copy_table_data(self):
        # Salin data dari tabel ke clipboard
        data = ""
        for row_id in self.table.get_children():
            row = self.table.item(row_id)['values']
            row_str = "\t".join(map(str, row))
            data += row_str + "\n"

        # Buat instance root untuk mengakses clipboard
        temp_root = tk.Tk()
        temp_root.withdraw()  # Sembunyikan jendela yang muncul
        temp_root.clipboard_clear()
        temp_root.clipboard_append(data)
        temp_root.update()  # Dibutuhkan untuk menyalin ke clipboard
        temp_root.destroy()  # Tutup root sementara

    def on_closing(self):
        #if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.destroy()

if __name__ == "__main__":
    app = HelixApp()
    app.mainloop()
