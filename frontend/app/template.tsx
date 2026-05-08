"use client";

import { motion } from "framer-motion";

// template.tsx re-renders on route change, so this gives every page a
// short cream-fade in. Keeps the rest of the layout (fonts, html shell)
// stable in layout.tsx so the transition feels like a clean wipe, not a
// full page reload.
export default function Template({ children }: { children: React.ReactNode }) {
  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.18, ease: [0.23, 1, 0.32, 1] }}
    >
      {children}
    </motion.div>
  );
}
