import type { Config } from "tailwindcss";

const config: Config = {
  content: [
    "./app/**/*.{ts,tsx}",
    "./components/**/*.{ts,tsx}",
    "./lib/**/*.{ts,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // GenomeScope warm palette - locked from project memory
        cream: {
          DEFAULT: "#F5F1E8",
          50: "#FBFAF5",
          100: "#F8F5EC",
          200: "#F5F1E8",
          300: "#EBE6D5",
          400: "#D9D2BC",
        },
        ink: {
          DEFAULT: "#1A1814",
          900: "#0F0E0B",
          800: "#1A1814",
          700: "#2D2A24",
          600: "#46423A",
          500: "#6B6657",
          400: "#928D7C",
          300: "#BEB9A7",
        },
        // Hero accent - yellow (Crextio-inspired)
        accent: {
          DEFAULT: "#F5D547",
          50: "#FEFCEC",
          100: "#FCF6C5",
          200: "#F9EC8A",
          300: "#F5D547",
          400: "#E8C220",
          500: "#C9A60E",
        },
        // Genomic feature colours - warm-family palette
        island: { DEFAULT: "#5BA89A", soft: "#D8EAE6" },     // muted teal
        hyper:  { DEFAULT: "#E89E3D", soft: "#F8E2C5" },     // amber
        promoter: { DEFAULT: "#7B3F61", soft: "#E5D2DD" },   // deep plum
      },
      fontFamily: {
        sans: ["var(--font-geist-sans)", "ui-sans-serif", "system-ui"],
        mono: ["var(--font-geist-mono)", "ui-monospace", "Menlo"],
      },
      borderRadius: {
        card: "20px",
        chip: "10px",
      },
      boxShadow: {
        card: "0 1px 2px rgba(26, 24, 20, 0.04), 0 8px 24px -12px rgba(26, 24, 20, 0.06)",
        cardHover:
          "0 1px 2px rgba(26, 24, 20, 0.04), 0 16px 40px -16px rgba(26, 24, 20, 0.10)",
      },
      transitionTimingFunction: {
        // Strong custom curves per emil-design-eng
        "out-strong": "cubic-bezier(0.23, 1, 0.32, 1)",
        "in-out-strong": "cubic-bezier(0.77, 0, 0.175, 1)",
      },
    },
  },
  plugins: [],
};
export default config;
